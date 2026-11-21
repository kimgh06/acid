"""1차 공정: 미강 → 효소 추출 → 고액분리 → Ca 침전 → 인산 전환 → 농축."""
import itertools
import math
from pathlib import Path

import warnings

import yaml
import biosteam as bst
import thermosteam as tmo
from biosteam.exceptions import CostWarning, DesignWarning

# 경제성은 범위 밖 → 장치 비용 상관식 경고 숨김
warnings.filterwarnings("ignore", category=CostWarning)
warnings.filterwarnings("ignore", category=DesignWarning)

import equilibrium as EQ
from species import P_FRACTION, P_MW, chemicals

PARAMS_FILE = Path(__file__).with_name("params.yaml")
SOLIDS = [c.ID for c in chemicals if c.locked_state == "s"]
SOLUTES = [c.ID for c in chemicals if c.ID not in SOLIDS and c.ID != "Water"]
P2O5_PER_H3PO4 = 141.9445 / (2 * 97.9952)
_run_ids = itertools.count()


class BalanceError(RuntimeError):
    """질량 또는 인 원소수지 불일치."""


# ---------------------------------------------------------------- 단위공정
class Tank(bst.Unit):
    """등온 반응/혼합조. 승온·반응열을 열유틸리티 부하로 계산."""
    _N_outs = 1
    _ins_size_is_fixed = False

    def _init(self, T, tau, reaction=None):
        self.T, self.tau, self.reaction = T, tau, reaction

    def _run(self):
        out = self.outs[0]
        out.mix_from(self.ins)
        if self.reaction:
            self.reaction(out)
        out.T = self.T

    def _design(self):
        duty = sum(s.Hnet for s in self.outs) - sum(s.Hnet for s in self.ins)
        self.add_heat_utility(duty, self.T)


class Separator(bst.Unit):
    """고액분리(여과·원심 공용). outs = (cake, liquid).

    cake 액분율만큼 액상이 같은 농도로 동반된다 → 용존 인 손실이 반영됨.
    ponytail: 평형 분리 모델. 세척·입도 영향은 2차.
    """
    _N_ins = 1
    _N_outs = 2

    def _init(self, recovery, moisture, bound=None):
        self.recovery, self.moisture, self.bound = recovery, moisture, bound or {}

    def _run(self):
        feed, (cake, liquid) = self.ins[0], self.outs
        cake.empty()
        cake.imol[SOLIDS] = feed.imol[SOLIDS] * self.recovery
        for ID, f in self.bound.items():  # 고상에 결합된 용존성분 (예: 미분해 피틴)
            cake.imol[ID] = feed.imol[ID] * f
        rest = feed.mol - cake.mol
        rest[[chemicals.index(i) for i in SOLIDS]] = 0  # 미회수 고형분은 액상 쪽으로
        solids_mass = cake.F_mass
        rest_mass = (rest * chemicals.MW).sum()
        need = solids_mass * self.moisture / (1 - self.moisture)
        cake.mol += rest * min(1.0, need / rest_mass if rest_mass else 0)
        liquid.mol = feed.mol - cake.mol
        cake.T = liquid.T = feed.T


class Washer(bst.Unit):
    """케이크 세척. ins = (cake, wash water), outs = (washed cake, wash liquor).

    케이크 함액 중 용질을 efficiency 만큼 회수하고, 빠진 자리는 물로 채운다(함액률 유지).
    ponytail: 다단 향류 세척을 회수율 하나로 요약. 단수·세척수 최적화는 2차.
    """
    _N_ins = 2
    _N_outs = 2

    def _init(self, efficiency):
        self.efficiency = efficiency

    def _run(self):
        (cake, water), (washed, liquor) = self.ins, self.outs
        washed.copy_like(cake)
        before = washed.F_mass
        washed.imol[SOLUTES] *= 1 - self.efficiency
        washed.imass["Water"] += before - washed.F_mass
        if washed.imass["Water"] > cake.imass["Water"] + water.imass["Water"] + 1e-9:
            raise ValueError(f"{self.ID}: 세척수 부족")
        liquor.mol = cake.mol + water.mol - washed.mol
        mix = tmo.Stream(None)
        mix.mix_from([cake, water])  # 단열 혼합온도 (세척수가 케이크를 식힘)
        washed.T = liquor.T = mix.T


# ---------------------------------------------------------------- 파라미터
FEEDSTOCKS = list(yaml.safe_load(PARAMS_FILE.read_text())["feedstocks"])


def load_params(overrides=None, feedstock=None):
    """params.yaml 로드. feedstock: 'feedstocks' 프리셋으로 feed 덮어쓰기.
    overrides = {'enzyme.T': 55, ...} → source='scenario'."""
    p = yaml.safe_load(PARAMS_FILE.read_text())
    presets = p.pop("feedstocks")
    if feedstock:
        p["feed"] = {**p["feed"], **presets[feedstock]}
    for key, value in (overrides or {}).items():
        group, name = key.split(".")
        p[group][name] = {**p[group][name], "value": value, "source": "scenario"}
    return p


def values(p):
    return {g: {k: v["value"] for k, v in d.items()} for g, d in p.items()}


def phytase_conversion(e, phytate_mol_per_kg):
    """피트산 → IP1 전환율. min(1차 동역학, 효소 처리용량). 식은 params.yaml 주석 참고.

    phytate_mol_per_kg: 미강 1 kg 당 피트산 mol.
    """
    f_T = math.exp(-(((e["T"] - e["T_opt"]) / e["T_sigma"]) ** 2))
    f_pH = math.exp(-(((e["pH"] - e["pH_opt"]) / e["pH_sigma"]) ** 2))
    kd = e["kd_slope"] * max(0.0, e["T"] - e["T_inact"])          # 열실활 1/h
    t_eff = e["time"] if kd == 0 else (1 - math.exp(-kd * e["time"])) / kd
    k = e["k_ref"] * (e["dose"] / e["dose_ref"]) * f_T * f_pH
    X_kinetic = 1 - math.exp(-k * t_eff)
    Pi_capacity = e["dose"] * 1e-6 * 60 * t_eff * e["activity_ratio"] * f_T * f_pH  # mol Pi/kg 미강
    X_capacity = Pi_capacity / (5 * phytate_mol_per_kg) if phytate_mol_per_kg else 1.0
    return {"X": min(X_kinetic, X_capacity, 1.0), "X_kinetic": X_kinetic, "X_capacity": X_capacity}


def feed_composition(f):
    """미강 kg/h 조성. 피트산·무기인은 총인과 유기인 비율에서 역산."""
    dry = f["bran_flow"] * (1 - f["moisture"])
    P = dry * f["total_P"]
    comp = {
        "Water": f["bran_flow"] - dry,
        "PhyticAcid": P * f["organic_P_frac"] / P_FRACTION["PhyticAcid"],
        "H3PO4": P * (1 - f["organic_P_frac"]) / P_FRACTION["H3PO4"],  # 무기인 대체
        "Protein": dry * f["protein"],
        "Lipid": dry * f["lipid"],
        "Carbohydrate": dry * f["carbohydrate"],
        "K_ion": dry * f["K"] * f["K_soluble"],      # 용출분만 용존, 나머지는 Ash 에 남음
        "Mg_ion": dry * f["Mg"] * f["Mg_soluble"],
    }
    comp["Ash"] = dry - sum(v for k, v in comp.items() if k != "Water")
    if comp["Ash"] < 0:
        raise ValueError(f"원료 조성 합 > 100 % (Ash = {comp['Ash']:.1f} kg/h)")
    return comp


# ---------------------------------------------------------------- 시스템
def build(params):
    v = values(params)
    bst.main_flowsheet.set_flowsheet(f"acid_{next(_run_ids)}")
    tmo.settings.set_thermo(chemicals)
    f, e = v["feed"], v["enzyme"]
    s1, pr, ac = v["separation1"], v["precipitation"], v["acidulation"]
    T = lambda c: c + 273.15

    bran = bst.Stream("bran", **feed_composition(f), units="kg/hr")
    water = bst.Stream("process_water", Water=f["bran_flow"] * v["pretreatment"]["water_ratio"], units="kg/hr")
    enzyme = bst.Stream("phytase", Phytase=f["bran_flow"] * e["dose"] / e["specific_activity"] / 1000, units="kg/hr")
    lime = bst.Stream("lime")
    acid = bst.Stream("sulfuric_acid")
    wash_water = bst.Stream("wash_water")

    enz = phytase_conversion(e, bran.imol["PhyticAcid"] * 1000 / f["bran_flow"])
    rx_phy = tmo.Reaction("PhyticAcid + 5 Water -> IP1 + 5 H3PO4", "PhyticAcid", enz["X"])
    ph = v["phosphatase"]
    X_ip1 = ph["conversion"] if ph["enabled"] else 0.0
    rx_ip1 = tmo.Reaction("IP1 + Water -> Inositol + H3PO4", "IP1", X_ip1)   # 2차 효소 (재순환 IP1 포함)
    rx_enz = tmo.SeriesReaction([rx_phy, rx_ip1])
    rx_prec = tmo.Reaction("2 H3PO4 + 3 CaOH2 -> TCP + 6 Water", "H3PO4", pr["conversion"])
    rx_acid_parts = [
        tmo.Reaction("TCP + 3 H2SO4 -> 2 H3PO4 + 3 Gypsum", "TCP", ac["conversion"]),
        tmo.Reaction("CaOH2 + H2SO4 -> Gypsum + 2 Water", "CaOH2", 1.0),  # 잔류 석회 중화
    ]
    rx_acid = tmo.ParallelReaction(rx_acid_parts)

    hr, r_rec = v["heat_recovery"], v["filtrate_recycle"]["fraction"]
    water_total = f["bran_flow"] * v["pretreatment"]["water_ratio"]
    # 여액 경로: S2 → [SP1 재순환 분기] → [HX1 공정수 예열] → 폐기(filtrate)
    #   r_rec > 0 이면 여액 일부를 R1 으로 되돌리고 신규 공정수는 부족분만 보충한다.
    filtrate_ID = "filtrate_hot" if (hr["enabled"] or r_rec) else "filtrate"
    HX1 = SP1 = MX0 = None
    water_to_R1 = water
    if hr["enabled"]:
        HX1 = bst.HXprocess("HX1_water_preheat", ins=(water, None), outs=("process_water_hot", "filtrate"),
                            dT=hr["dT_min"], phase0="l", phase1="l")
        water_to_R1 = HX1 - 0
    if r_rec:
        SP1 = bst.Splitter("SP1_filtrate_split", ins=None, outs=("filtrate_recycle", "filtrate_purge" if HX1 else "filtrate"),
                           split=r_rec)
        MX0 = bst.Mixer("MX0_water_mix", ins=(water_to_R1, SP1 - 0), outs="water_to_R1")
        water_to_R1 = MX0 - 0
    R1 = Tank("R1_extraction", ins=(bran, water_to_R1, enzyme), T=T(e["T"]), tau=e["time"], reaction=rx_enz)
    S1 = Separator("S1_solid_liquid", ins=R1 - 0, outs=("residue_cake", "extract"),
                   recovery=s1["solids_recovery"], moisture=s1["cake_moisture"],
                   bound={"PhyticAcid": s1["phytate_to_cake"]})
    R2 = Tank("R2_precipitation", ins=(S1 - 1, lime), T=T(pr["T"]), tau=pr["time"], reaction=rx_prec)
    S2 = Separator("S2_TCP_filter", ins=R2 - 0, outs=("TCP_cake", filtrate_ID),
                   recovery=1.0, moisture=pr["cake_moisture"])
    R3 = Tank("R3_acidulation", ins=(S2 - 0, acid), T=T(ac["T"]), tau=ac["time"], reaction=rx_acid)
    S3 = Separator("S3_gypsum_filter", ins=R3 - 0, outs=("gypsum_wet", "crude_acid"),
                   recovery=1.0, moisture=ac["cake_moisture"])
    W1 = Washer("W1_gypsum_wash", ins=(S3 - 0, wash_water), outs=("gypsum_cake", "wash_liquor"),
                efficiency=ac["wash_efficiency"])
    M1 = bst.Mixer("M1_acid_mix", ins=(S3 - 1, W1 - 1), outs="acid_to_E1")  # ponytail: 세척액을 반응기로 재순환하지 않고 농축으로 보냄
    E1 = bst.MultiEffectEvaporator("E1_concentration", ins=M1 - 0, outs=("product_acid", "condensate"),
                                   P=(v["concentration"]["pressure"],), V=0.5, V_definition="Overall")

    @R2.add_specification(run=True)
    def dose_lime():
        ext = S1.outs[1]
        lime.imol["CaOH2"] = pr["Ca_P_ratio"] * ext.imol["H3PO4"]
        # Mg 공침분: 추출액 농도에서 PHREEQC 평형 (pH ≥ 9 에서 브루사이트) → S2 케이크로
        kgw = ext.imass["Water"] / 1000
        if kgw <= 0 or ext.imol["Mg_ion"] <= 0 or ext.imol["H3PO4"] <= 0:  # 재순환 첫 반복 등 빈 추출액
            S2.bound["Mg_ion"] = 0.0
            return
        P_mM, Mg_mM = ext.imol["H3PO4"] / kgw * 1000, ext.imol["Mg_ion"] / kgw * 1000
        S2.bound["Mg_ion"] = EQ.removal(P_mM, pr["Ca_P_ratio"] * P_mM, pr["pH"], pr["T"], "TCP", Mg_mM)["Mg_removal"]

    @R3.add_specification(run=True)
    def dose_acid():
        cake = S2.outs[0]
        acid.imol["H2SO4"] = ac["acid_excess"] * (3 * cake.imol["TCP"] + cake.imol["CaOH2"])
        acid.imass["Water"] = acid.imass["H2SO4"] * (1 - ac["acid_strength"]) / ac["acid_strength"]

    @W1.add_specification(run=True)
    def dose_wash_water():
        cake = S3.outs[0]
        wash_water.imass["Water"] = ac["wash_water"] * (cake.F_mass - cake.imass[SOLIDS].sum())

    @E1.add_specification(run=True)
    def target_concentration():
        feed = E1.ins[0]
        w = v["concentration"]["target_wt"]
        water_keep = feed.imass["H3PO4"] / w - (feed.F_mass - feed.imass["Water"])
        if water_keep < 0:
            raise ValueError(f"불순물이 많아 목표 농도 {w:.0%} 도달 불가")
        if water_keep >= feed.imass["Water"]:
            raise ValueError(f"증발기 입구가 이미 목표 농도 {w:.0%} 이상 — 농축 불필요")
        # V: VLE 대상(=물, 나머지는 상 고정) 몰 기준 증발분율
        E1.V = 1 - water_keep / feed.imass["Water"]

    units = [R1, S1, R2, S2, R3, S3, W1, M1, E1]
    if SP1:
        SP1.ins[0] = S2.outs[1]
        units[:0] = [MX0, SP1]

        def makeup_water():  # 신규 공정수 = 목표 물량 - 재순환 물
            fresh = water_total - SP1.outs[0].imass["Water"]
            if fresh < 0:
                raise ValueError("재순환 물이 목표 물량 초과 — filtrate_recycle.fraction 을 낮춘다")
            water.imass["Water"] = fresh
        (HX1 or MX0).add_specification(makeup_water, run=True)
    if HX1:
        HX1.ins[1] = SP1.outs[1] if SP1 else S2.outs[1]
        units.insert(0, HX1)
    sys = bst.main_flowsheet.create_system("acid_process")  # 재순환 자동 인식
    # 재순환 수렴을 수지 검사(상대 1e-6)보다 엄격하게. 기본값(~1e-3)이면 R1 에서 BalanceError.
    sys.set_tolerance(mol=1e-9, rmol=1e-10, T=1e-6, rT=1e-9, subsystems=True, maxiter=500)
    enz["X_IP1"] = X_ip1
    return sys, units, {"enzyme": enz, "reactions": (rx_phy, rx_ip1, rx_prec, *rx_acid_parts)}


# ---------------------------------------------------------------- 수지·결과
def P_flow(s):
    return s.get_atomic_flow("P") * P_MW  # kg P/h


ENERGY_TOL_KJ_H = 1.0  # 절대 허용오차. Hnet(생성열 포함)이 ~1e8 이라 상대오차는 오류를 가린다.


def energy_residual(u):
    """입열 + 유틸리티 실부하(unit_duty, 전달효율 손실 제외) - 출열, kJ/h."""
    q = sum(hu.unit_duty for hu in u.heat_utilities)
    return sum(s.Hnet for s in u.ins) + q - sum(s.Hnet for s in u.outs)


def check_balance(units, tol=1e-6):
    """단위공정별 질량·P 원소수지(상대오차 > tol)와 에너지 수지(|잔차| > 1 kJ/h) → BalanceError."""
    for u in units:
        for name, f in (("mass", lambda s: s.F_mass), ("P", P_flow)):
            a, b = sum(f(s) for s in u.ins), sum(f(s) for s in u.outs)
            if abs(a - b) > tol * max(abs(a), 1e-12):
                raise BalanceError(f"{u.ID} {name} 수지 불일치: in {a:.6g} / out {b:.6g}")
        e = energy_residual(u)
        if abs(e) > ENERGY_TOL_KJ_H:
            raise BalanceError(f"{u.ID} energy 수지 불일치: 잔차 {e:.3g} kJ/h")


def run(overrides=None, feedstock=None):
    params = load_params(overrides, feedstock)
    sys, units, info = build(params)
    for rx in info["reactions"]:
        rx.check_atomic_balance()
    sys.simulate()
    check_balance(units)

    F = bst.main_flowsheet.stream
    feeds = [F.bran, F.process_water, F.phytase, F.lime, F.sulfuric_acid, F.wash_water]
    product, E1 = F.product_acid, units[-1]
    P_in = P_flow(F.bran)
    dry = product.F_mass - product.imass["Water"]
    heat = [hu.duty for u in units for hu in u.heat_utilities]
    mass_in, mass_out = sum(s.F_mass for s in feeds), sum(s.F_mass for s in sys.products)
    steam_kg_h = sum(hu.flow * 18.01528 for hu in E1.heat_utilities if hu.duty > 0)
    enz = info["enzyme"]

    result = {
        "phytate_conversion": enz["X"],                    # 피트산 → IP1
        "phytate_P_released": enz["X"] * (5 + enz["X_IP1"]) / 6,  # 피트산 P 중 무기인으로 방출 (2차 효소 포함, 1회 통과)
        "enzyme_limited_by": "capacity" if enz["X_capacity"] < enz["X_kinetic"] else "kinetics",
        "P_extraction": 1 - P_flow(F.residue_cake) / P_in,  # 순 추출률: 잔사로 빠지지 않은 원료 P (재순환 무관)
        "P_recovery": P_flow(product) / P_in,              # 원료 대비 총인 회수율 (제품 중 모든 P)
        "H3PO4_P_recovery": product.imass["H3PO4"] * P_FRACTION["H3PO4"] / P_in,  # H3PO4 형태 P 만
        "H3PO4_kg_h": product.imass["H3PO4"],
        "product_kg_h": product.F_mass,
        "H3PO4_wt": product.imass["H3PO4"] / product.F_mass,
        "purity_dry": product.imass["H3PO4"] / dry,        # 무수 기준 순도
        **{f"P_loss_{s.ID}_kg_h": P_flow(s) for s in (F.residue_cake, F.filtrate, F.gypsum_cake, F.condensate)},
        "acidulation_P_recovery": P_flow(F.acid_to_E1) / P_flow(F.TCP_cake),  # 전환+세척 후 산으로 회수
        "water_kg_h": F.process_water.F_mass + F.sulfuric_acid.imass["Water"] + F.wash_water.F_mass,
        "lime_kg_h": F.lime.F_mass,
        "H2SO4_kg_h": F.sulfuric_acid.imass["H2SO4"],
        "phytase_kg_h": F.phytase.F_mass,
        "heating_MJ_h": sum(d for d in heat if d > 0) / 1000,
        "cooling_MJ_h": -sum(d for d in heat if d < 0) / 1000,
        "electricity_kW": sum(u.power_utility.rate for u in units),
        "E1_steam_t_per_t_P2O5": steam_kg_h / (product.imass["H3PO4"] * P2O5_PER_H3PO4),
        "E1_steam_t_per_t_water": steam_kg_h / F.condensate.F_mass,
        "E1_feed_P2O5_wt": F.acid_to_E1.imass["H3PO4"] * P2O5_PER_H3PO4 / F.acid_to_E1.F_mass,
        "mass_balance_error": (mass_out - mass_in) / mass_in,
        "energy_residual_max_kJ_h": max(abs(energy_residual(u)) for u in units),
        "gypsum_kg_h": F.gypsum_cake.imass["Gypsum"],
        "IP1_kg_h": sum(s.imass["IP1"] for s in sys.products),
        "inositol_kg_h": sum(s.imass["Inositol"] for s in sys.products),
        "Mg_to_TCP_cake": units[[u.ID for u in units].index("S2_TCP_filter")].bound.get("Mg_ion", 0.0),
        "product_K_wt": product.imass["K_ion"] / product.F_mass,
        "product_Mg_wt": product.imass["Mg_ion"] / product.F_mass,
    }
    return result, sys, params


def utility_table(sys, result):
    """유틸리티·약품 사용량과 P2O5 1 t 당 원단위 (스팀·냉각수는 에이전트 유량, 전달효율 손실 포함)."""
    import pandas as pd
    F = sys.flowsheet.stream
    P2O5 = F.product_acid.imass["H3PO4"] * P2O5_PER_H3PO4
    hus = [hu for u in sys.units for hu in u.heat_utilities]
    rows = {
        "steam_kg_h": sum(hu.flow * 18.01528 for hu in hus if hu.duty > 0),
        "cooling_water_kg_h": sum(hu.flow * 18.01528 for hu in hus if hu.duty < 0),
        "power_kW": result["electricity_kW"],
        "water_kg_h": result["water_kg_h"],
        "lime_kg_h": result["lime_kg_h"],
        "H2SO4_kg_h": result["H2SO4_kg_h"],
        "phytase_kg_h": result["phytase_kg_h"],
    }
    # 원단위: 질량 유량은 t/t P2O5, 전력은 kWh/t P2O5
    return pd.DataFrame([{"item": k, "use": v, "per_t_P2O5": v / P2O5 * (1000 if k == "power_kW" else 1),
                          "per_t_unit": "kWh/t" if k == "power_kW" else "t/t"} for k, v in rows.items()]
                        + [{"item": "P2O5_kg_h", "use": P2O5, "per_t_P2O5": 1.0, "per_t_unit": "t/t"}]).set_index("item")


def stream_table(sys):
    import pandas as pd
    rows = [{
        "stream": s.ID,
        "from": s.source.ID if s.source else "-",
        "to": s.sink.ID if s.sink else "-",
        "kg_h": s.F_mass,
        "P_kg_h": P_flow(s),
        "H3PO4_kg_h": s.imass["H3PO4"],
        "T_C": s.T - 273.15,
    } for s in sys.streams]
    return pd.DataFrame(rows).set_index("stream")


def duty_table(units):
    import pandas as pd
    return pd.DataFrame([{
        "unit": u.ID,
        "heating_MJ_h": sum(hu.duty for hu in u.heat_utilities if hu.duty > 0) / 1000,
        "cooling_MJ_h": -sum(hu.duty for hu in u.heat_utilities if hu.duty < 0) / 1000,
        "power_kW": u.power_utility.rate,
    } for u in units]).set_index("unit")


def values_by_source(params, source):
    """source('assumed' / 'literature' / 'measured' / 'scenario') 인 항목 목록."""
    return [f"{g}.{k} = {v['value']} {v['unit']}" + (f"  [{v['ref']}]" if v.get("ref") else "")
            for g, d in params.items() for k, v in d.items() if v["source"] == source]


def assumed_values(params):
    return values_by_source(params, "assumed")
