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

from species import P_FRACTION, P_MW, chemicals

PARAMS_FILE = Path(__file__).with_name("params.yaml")
SOLIDS = [c.ID for c in chemicals if c.locked_state == "s"]
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


# ---------------------------------------------------------------- 파라미터
def load_params(overrides=None):
    """params.yaml 로드. overrides = {'enzyme.T': 55, ...} → source='scenario'."""
    p = yaml.safe_load(PARAMS_FILE.read_text())
    for key, value in (overrides or {}).items():
        group, name = key.split(".")
        p[group][name] = {**p[group][name], "value": value, "source": "scenario"}
    return p


def values(p):
    return {g: {k: v["value"] for k, v in d.items()} for g, d in p.items()}


def phytase_conversion(e):
    """X = 1 - exp(-k t). 가정 모델 — 문헌·실험 동역학으로 교체 대상."""
    f_T = math.exp(-(((e["T"] - e["T_opt"]) / e["T_sigma"]) ** 2))
    f_pH = math.exp(-(((e["pH"] - e["pH_opt"]) / e["pH_sigma"]) ** 2))
    k = e["k_ref"] * (e["dose"] / e["dose_ref"]) * f_T * f_pH
    return 1 - math.exp(-k * e["time"])


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

    X_enz = phytase_conversion(e)
    rx_enz = tmo.Reaction("PhyticAcid + 6 Water -> Inositol + 6 H3PO4", "PhyticAcid", X_enz)
    rx_prec = tmo.Reaction("2 H3PO4 + 3 CaOH2 -> TCP + 6 Water", "H3PO4", pr["conversion"])
    rx_acid_parts = [
        tmo.Reaction("TCP + 3 H2SO4 -> 2 H3PO4 + 3 Gypsum", "TCP", ac["conversion"]),
        tmo.Reaction("CaOH2 + H2SO4 -> Gypsum + 2 Water", "CaOH2", 1.0),  # 잔류 석회 중화
    ]
    rx_acid = tmo.ParallelReaction(rx_acid_parts)

    R1 = Tank("R1_extraction", ins=(bran, water, enzyme), T=T(e["T"]), tau=e["time"], reaction=rx_enz)
    S1 = Separator("S1_solid_liquid", ins=R1 - 0, outs=("residue_cake", "extract"),
                   recovery=s1["solids_recovery"], moisture=s1["cake_moisture"],
                   bound={"PhyticAcid": s1["phytate_to_cake"]})
    R2 = Tank("R2_precipitation", ins=(S1 - 1, lime), T=T(pr["T"]), tau=pr["time"], reaction=rx_prec)
    S2 = Separator("S2_TCP_filter", ins=R2 - 0, outs=("TCP_cake", "filtrate"),
                   recovery=1.0, moisture=pr["cake_moisture"])
    R3 = Tank("R3_acidulation", ins=(S2 - 0, acid), T=T(ac["T"]), tau=ac["time"], reaction=rx_acid)
    S3 = Separator("S3_gypsum_filter", ins=R3 - 0, outs=("gypsum_cake", "crude_acid"),
                   recovery=1.0, moisture=ac["cake_moisture"])
    E1 = bst.MultiEffectEvaporator("E1_concentration", ins=S3 - 1, outs=("product_acid", "condensate"),
                                   P=(101325, 60000, 25000), V=0.5, V_definition="Overall")

    @R2.add_specification(run=True)
    def dose_lime():
        lime.imol["CaOH2"] = 1.5 * pr["lime_excess"] * S1.outs[1].imol["H3PO4"]

    @R3.add_specification(run=True)
    def dose_acid():
        cake = S2.outs[0]
        acid.imol["H2SO4"] = ac["acid_excess"] * (3 * cake.imol["TCP"] + cake.imol["CaOH2"])
        acid.imass["Water"] = acid.imass["H2SO4"] * (1 - ac["acid_strength"]) / ac["acid_strength"]

    @E1.add_specification(run=True)
    def target_concentration():
        feed = E1.ins[0]
        w = v["concentration"]["target_wt"]
        water_keep = feed.imass["H3PO4"] / w - (feed.F_mass - feed.imass["Water"])
        if water_keep < 0:
            raise ValueError(f"불순물이 많아 목표 농도 {w:.0%} 도달 불가")
        # V: VLE 대상(=물, 나머지는 상 고정) 몰 기준 증발분율
        E1.V = 1 - water_keep / feed.imass["Water"]

    units = [R1, S1, R2, S2, R3, S3, E1]
    sys = bst.System("acid_process", path=units)
    return sys, units, {"X_enzyme": X_enz, "reactions": (rx_enz, rx_prec, *rx_acid_parts)}


# ---------------------------------------------------------------- 수지·결과
def P_flow(s):
    return s.get_atomic_flow("P") * P_MW  # kg P/h


def check_balance(units, tol=1e-6):
    """단위공정별 전체 질량·인 원소수지. 상대오차 > tol 이면 BalanceError."""
    for u in units:
        for name, f in (("mass", lambda s: s.F_mass), ("P", P_flow)):
            a, b = sum(f(s) for s in u.ins), sum(f(s) for s in u.outs)
            if abs(a - b) > tol * max(abs(a), 1e-12):
                raise BalanceError(f"{u.ID} {name} 수지 불일치: in {a:.6g} / out {b:.6g}")


def run(overrides=None):
    params = load_params(overrides)
    sys, units, info = build(params)
    for rx in info["reactions"]:
        rx.check_atomic_balance()
    sys.simulate()
    check_balance(units)

    F = bst.main_flowsheet.stream
    feeds, product = [F.bran, F.process_water, F.phytase, F.lime, F.sulfuric_acid], F.product_acid
    P_in = P_flow(F.bran)
    dry = product.F_mass - product.imass["Water"]
    heat = [hu.duty for u in units for hu in u.heat_utilities]
    mass_in, mass_out = sum(s.F_mass for s in feeds), sum(s.F_mass for s in sys.products)

    result = {
        "phytate_conversion": info["X_enzyme"],
        "P_extraction": P_flow(F.extract) / P_in,          # 액상으로 추출된 인
        "P_recovery": P_flow(product) / P_in,              # 원료 대비 총인 회수율
        "H3PO4_kg_h": product.imass["H3PO4"],
        "product_kg_h": product.F_mass,
        "H3PO4_wt": product.imass["H3PO4"] / product.F_mass,
        "purity_dry": product.imass["H3PO4"] / dry,        # 무수 기준 순도
        **{f"P_loss_{s.ID}_kg_h": P_flow(s) for s in (F.residue_cake, F.filtrate, F.gypsum_cake, F.condensate)},
        "water_kg_h": F.process_water.F_mass + F.sulfuric_acid.imass["Water"],
        "lime_kg_h": F.lime.F_mass,
        "H2SO4_kg_h": F.sulfuric_acid.imass["H2SO4"],
        "phytase_kg_h": F.phytase.F_mass,
        "heating_MJ_h": sum(d for d in heat if d > 0) / 1000,
        "cooling_MJ_h": -sum(d for d in heat if d < 0) / 1000,
        "electricity_kW": sum(u.power_utility.rate for u in units),
        "mass_balance_error": (mass_out - mass_in) / mass_in,
        "gypsum_kg_h": F.gypsum_cake.imass["Gypsum"],
        "inositol_kg_h": sum(s.imass["Inositol"] for s in sys.products),
    }
    return result, sys, params


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


def assumed_values(params):
    return [f"{g}.{k} = {v['value']} {v['unit']}" for g, d in params.items()
            for k, v in d.items() if v["source"] == "assumed"]
