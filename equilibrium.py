"""PHREEQC 침전 평형. VAL-03 물성 검증 + R2 Mg 공침 분율 계산 (process.dose_lime).

ponytail: pH 는 입력(Fix_H+ 로 고정, NaOH 로 조절). 석회 소요량과 pH 의 연동은 미계산.
"""
import re
from pathlib import Path

import phreeqc
from phreeqc import Phreeqc

DB = "minteq.v4.dat"
MINERALS = {"HAP": "Hydroxylapatite", "TCP": "Ca3(PO4)2(beta)", "brushite": "CaHPO4:2H2O"}
MG_MINERALS = ("Brucite", "Mg3(PO4)2", "MgHPO4:3H2O")


def _log_k(reaction):
    """minteq.v4.dat 에서 반응식 바로 다음 줄의 log_k."""
    text = (Path(phreeqc.__file__).parent / "databases" / DB).read_text(errors="ignore")
    m = re.search(re.escape(reaction) + r"\s*\n\s*log_k\s+(-?[\d.]+)", text)
    return float(m.group(1))


def db_constants():
    """DB 값에서 인산 pKa 와 HAP log Ksp (Ca5(PO4)3OH = 5Ca + 3PO4 + OH 형식) 계산."""
    k1, k2, k3 = (_log_k(r) for r in ("H+ + PO4-3 = HPO4-2", "2 H+ + PO4-3 = H2PO4-", "3 H+ + PO4-3 = H3PO4"))
    kw = _log_k("H2O = OH- + H+")
    hap = _log_k("Ca5(PO4)3OH + H+ = 5 Ca+2 + 3 PO4-3 + H2O")
    return {"pKa1": k3 - k2, "pKa2": k2 - k1, "pKa3": k1, "HAP_logKsp": hap + kw}


def removal(P_mM, Ca_mM, pH, T_C, mineral="HAP", Mg_mM=0.0):
    """Ca-인산 mineral (+ Mg_mM > 0 이면 Mg 광물 3종) 평형 침전 시 P·Mg 제거율 (0~1)."""
    mg_lines = "".join(f"\n    {m} 0 0" for m in MG_MINERALS) if Mg_mM else ""
    p = Phreeqc()
    p.LoadBuiltInDatabase(DB)
    p.RunString(f"""
PHASES
Fix_H+
    H+ = H+
    log_k 0
SOLUTION 1
    units mmol/kgw
    temp {T_C}
    pH 7
    P {P_mM}
    Ca {Ca_mM}
    Mg {Mg_mM}
    Cl 1 charge
EQUILIBRIUM_PHASES 1
    Fix_H+ {-pH} NaOH 10
    {MINERALS[mineral]} 0 0{mg_lines}
SELECTED_OUTPUT
    -reset false
    -pH true
    -totals P Mg
END
""")
    out = p.GetSelectedOutput()
    P_left = out["P(mol/kgw)"][-1] * 1000
    res = {"removal": 1 - P_left / P_mM, "pH": out["pH"][-1]}
    if Mg_mM:
        res["Mg_removal"] = 1 - out["Mg(mol/kgw)"][-1] * 1000 / Mg_mM
    return res


if __name__ == "__main__":
    c = db_constants()
    print({k: round(v, 3) for k, v in c.items()})
    for m in MINERALS:
        print(m, [round(removal(55, 110, pH, 50, m)["removal"], 4) for pH in (5, 6, 7, 8, 8.5, 9, 10)])
