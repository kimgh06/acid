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


def test_imbalance_is_detected():
    _, sys, _ = P.run()
    sys.flowsheet.stream.extract.imass["H3PO4"] += 1  # 인 1 kg/h 증발시킴
    try:
        P.check_balance(sys.units)
    except P.BalanceError:
        return
    raise AssertionError("수지 불일치를 못 잡음")


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


def test_scenario_marked():
    _, _, params = P.run({"enzyme.T": 55})
    assert params["enzyme"]["T"]["source"] == "scenario"
    assert "enzyme.T = 55" not in " ".join(P.assumed_values(params))


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
            print("ok", name)
