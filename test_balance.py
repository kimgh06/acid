"""VAL-01 수지 검증 + VAL-04 조건 변화 타당성. 실행: uv run test_balance.py (pytest도 가능)"""
import process as P


def test_balance_closes():
    r, sys, _ = P.run()
    assert abs(r["mass_balance_error"]) < 1e-9
    F = sys.flowsheet.stream
    P_in = P.P_flow(F.bran)
    P_out = sum(P.P_flow(s) for s in sys.products)
    assert abs(P_out - P_in) < 1e-9 * P_in
    # 회수 + 공정별 손실 = 투입
    losses = sum(v for k, v in r.items() if k.startswith("P_loss_"))
    assert abs(r["P_recovery"] * P_in + losses - P_in) < 1e-9 * P_in
    assert 0 < r["P_recovery"] < r["P_extraction"] <= 1
    assert r["H3PO4_P_recovery"] <= r["P_recovery"]     # 제품 중 유기 P(IP1) 제외


def test_imbalance_is_detected():
    _, sys, _ = P.run()
    sys.flowsheet.stream.extract.imass["H3PO4"] += 1  # 인 1 kg/h 증발시킴
    try:
        P.check_balance(sys.units)
    except P.BalanceError:
        return
    raise AssertionError("수지 불일치를 못 잡음")


def test_energy_balance():
    r, sys, _ = P.run()
    assert r["energy_residual_max_kJ_h"] < P.ENERGY_TOL_KJ_H
    sys.flowsheet.stream.extract.T += 5  # 추출액만 5 K 올림 → S1 에너지 불일치
    try:
        P.check_balance(sys.units)
    except P.BalanceError as e:
        assert "energy" in str(e)
        return
    raise AssertionError("에너지 수지 불일치를 못 잡음")


def test_utility_table():
    r, sys, _ = P.run()
    u = P.utility_table(sys, r)
    P2O5 = u.loc["P2O5_kg_h", "use"]
    assert abs(u.loc["power_kW", "per_t_P2O5"] - r["electricity_kW"] / P2O5 * 1000) < 1e-9   # kWh/t
    assert abs(u.loc["lime_kg_h", "per_t_P2O5"] - r["lime_kg_h"] / P2O5) < 1e-9


def test_target_concentration():
    for w in (0.6, 0.85):
        r, *_ = P.run({"concentration.target_wt": w})
        assert abs(r["H3PO4_wt"] - w) < 1e-6


def test_condition_trends():
    # dose 500: 효소용량이 제한하는 영역에서 경향 확인 (기준 2500은 거의 100% 포화)
    X = lambda **kw: P.run({"enzyme.dose": 500, **{f"enzyme.{k}": v for k, v in kw.items()}})[0]["phytate_conversion"]
    assert X(time=1) < X(time=2) < X(time=4)                # 시간 ↑ → 전환율 ↑
    assert X(pH=5.5) > X(pH=3.5) and X(pH=5.5) > X(pH=7.5)  # 최적 pH
    assert X(dose=1000) > X()                               # 투입량 ↑
    assert X(T=70) < X(T=55)                                # 55°C 초과 열실활


def test_phytase_endpoint():
    r, *_ = P.run({"enzyme.dose": 1e6})                     # 효소 과잉이어도
    assert r["phytate_conversion"] <= 1
    assert abs(r["phytate_P_released"] - 5 / 6) < 1e-3     # IP1에서 멈춤 → 최대 5/6


def test_wash_reduces_loss():
    washed, *_ = P.run()
    unwashed, *_ = P.run({"acidulation.wash_efficiency": 0.0})
    assert washed["P_loss_gypsum_cake_kg_h"] < unwashed["P_loss_gypsum_cake_kg_h"]
    assert washed["P_recovery"] > unwashed["P_recovery"]


def test_phreeqc_constants():
    import equilibrium as EQ
    c = EQ.db_constants()
    for k, ref in (("pKa1", 2.16), ("pKa2", 7.21), ("pKa3", 12.32)):
        assert abs(c[k] - ref) < 0.1, (k, c[k])
    assert abs(c["HAP_logKsp"] - (-58.33)) < 0.2         # McDowell 1977
    lo, hi = (EQ.removal(55, 110, pH, 50, "TCP")["removal"] for pH in (5, 8))
    assert lo < hi                                       # pH ↑ → 침전 ↑


def test_feedstocks():
    for fs in P.FEEDSTOCKS:
        r, _, params = P.run(feedstock=fs)                # 수지 검사 포함
        assert 0 < r["P_recovery"] < 1
        assert params["feed"]["lipid"]["value"] != P.load_params()["feed"]["lipid"]["value"]


def test_uncertainty_ranges():
    import uncertainty as U
    params = P.load_params()
    for key, (lo, hi) in U.ranged(params).items():
        g, k = key.split(".")
        assert lo <= params[g][k]["value"] <= hi, key       # 기준값이 범위 안
    df, _, _ = U.sample(3)
    assert len(df) == 3


def test_filtrate_recycle():
    base, *_ = P.run()
    for hr in (1, 0):                                   # 열회수 유무 두 배관 모두
        rec, sys, _ = P.run({"filtrate_recycle.fraction": 0.5, "heat_recovery.enabled": hr})   # 수지 검사 포함
        assert rec["water_kg_h"] < base["water_kg_h"]
        assert 0 < rec["P_recovery"] < rec["P_extraction"] <= 1   # 재순환 P 를 추출로 이중 계산하지 않음
        F = sys.flowsheet.stream
        assert abs(F.filtrate.F_mass - 0.5 * F.filtrate_hot.F_mass) < 1e-6 * F.filtrate_hot.F_mass  # 방류 = (1-r) × 여액


def test_mg_k_tracking():
    base, sys, _ = P.run()
    F = sys.flowsheet.stream
    for ion in ("K_ion", "Mg_ion"):                     # 금속 보존: 원료 용출분 = 배출 합
        assert abs(F.bran.imass[ion] - sum(s.imass[ion] for s in sys.products)) < 1e-9 * F.bran.imass[ion]
    high, *_ = P.run({"precipitation.pH": 9.0})
    assert base["Mg_to_TCP_cake"] < 0.05 < 0.5 < high["Mg_to_TCP_cake"]   # pH 9 브루사이트 공침
    assert high["product_Mg_wt"] > 10 * base["product_Mg_wt"]
    assert high["purity_dry"] < base["purity_dry"]


def test_phosphatase():
    base, *_ = P.run()
    zero, *_ = P.run({"phosphatase.enabled": 1, "phosphatase.conversion": 0.0})
    on, sys, _ = P.run({"phosphatase.enabled": 1})                     # 원소수지·수지 검사 포함
    assert abs(zero["P_recovery"] - base["P_recovery"]) < 1e-9        # 전환 0 = 미사용
    assert on["P_recovery"] > base["P_recovery"] + 0.05
    assert on["IP1_kg_h"] < base["IP1_kg_h"] and on["inositol_kg_h"] > 0
    assert on["phytate_P_released"] > 5 / 6                           # IP1 종점 상한 돌파


def test_scenario_marked():
    _, _, params = P.run({"enzyme.T": 55})
    assert params["enzyme"]["T"]["source"] == "scenario"
    assert "enzyme.T = 55" not in " ".join(P.assumed_values(params))


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
            print("ok", name)
