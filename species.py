"""공정 물질 정의. DB에 없는 원료 성분은 대체물질(surrogate)로 둔다."""
import warnings

warnings.filterwarnings("ignore", category=SyntaxWarning)  # biosteam docstring 수식 경고

import thermosteam as tmo

P_MW = 30.973762


def _lock(c, phase):
    if c.locked_state != phase:
        c.at_state(phase)  # 제자리 고정, None 반환
    return c


def _chem(ID, search=None, phase=None):
    c = tmo.Chemical(ID, search_ID=search) if search else tmo.Chemical(ID)
    return _lock(c, phase) if phase else c


def _surrogate(ID, base, phase):
    return _lock(tmo.Chemical(base).copy(ID), phase)


# 수용성 비휘발 성분은 'l'로 고정해 증발기에서 날아가지 않게 한다.
# ponytail: 이온 평형 없음 — 해리·pH는 PHREEQC 연동 시 다룬다.
chemicals = tmo.Chemicals([
    "Water",
    _chem("PhyticAcid", phase="l"),           # 피트산(피틴 대체), C6H18O24P6
    _chem("Inositol", "myo-inositol", phase="l"),
    _chem("H3PO4", phase="l"),
    _chem("H2SO4", phase="l"),
    _chem("CaOH2", "Ca(OH)2", phase="s"),
    _chem("TCP", "Ca3(PO4)2", phase="s"),     # 인산삼칼슘
    _chem("Gypsum", "CaSO4", phase="s"),      # ponytail: 무수물 — 실제 습식법은 CaSO4·2H2O
    # 대체물질: 인 수지에 영향 없음, 질량 운반용
    _surrogate("Protein", "Cellulose", "s"),
    _surrogate("Lipid", "Triolein", "s"),
    _surrogate("Carbohydrate", "Cellulose", "s"),
    _surrogate("Ash", "SiO2", "s"),           # 인 외 무기질
    _surrogate("Phytase", "Glucose", "l"),  # 효소 단백질(수용성)
])

# DB의 피트산 Hf(+9004 kJ/mol)는 기여법 오류값. 가수분해 반응열 = 0 이 되도록 맞춘다 (가정값).
_h = {c.ID: c.Hf for c in chemicals}
chemicals.PhyticAcid.Hf = _h["Inositol"] + 6 * _h["H3PO4"] - 6 * _h["Water"]
for ID in ("Protein", "Lipid", "Carbohydrate", "Ash", "Phytase"):
    chemicals[ID].Hf = 0  # 비반응 대체물질

chemicals.compile()
chemicals.set_synonym("Water", "H2O")
tmo.settings.set_thermo(chemicals)

# 원소 P 질량분율
P_FRACTION = {c.ID: c.atoms.get("P", 0) * P_MW / c.MW for c in chemicals}
