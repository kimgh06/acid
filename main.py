"""기준 공정 결과 + 조건별 비교(SIM-09, SIM-10) + 문헌·물성 검증(VAL-02, VAL-03).

실행: uv run main.py [--out DIR]   (--out: 모든 표를 DIR 에 CSV·JSON 으로 저장, 보고서 부록용)
"""
import json
import sys as _sys
from pathlib import Path

import pandas as pd

import equilibrium as EQ
import process as P

SCENARIOS = {
    "base": {},
    "dose 1000": {"enzyme.dose": 1000},
    "pH 4.5": {"enzyme.pH": 4.5},
    "pH 6.5": {"enzyme.pH": 6.5},
    "T 45C": {"enzyme.T": 45},
    "T 65C": {"enzyme.T": 65},          # 열실활 영역
    "time 1h": {"enzyme.time": 1.0},
    "Ca/P 1.67": {"precipitation.Ca_P_ratio": 1.67},
    "no wash": {"acidulation.wash_efficiency": 0.0},
    "TCP cake 30%": {"precipitation.cake_moisture": 0.3},   # 필터프레스 탈수 가정
    "no heat recovery": {"heat_recovery.enabled": 0},
    "recycle 50%": {"filtrate_recycle.fraction": 0.5},
    "recycle 90%": {"filtrate_recycle.fraction": 0.9},
    "precip pH 9": {"precipitation.pH": 9.0},          # 브루사이트 공침 → 제품 Mg
    "phosphatase 50%": {"phosphatase.enabled": 1, "phosphatase.conversion": 0.5},
    "phosphatase 90%": {"phosphatase.enabled": 1},
    "phosphatase+recycle 50%": {"phosphatase.enabled": 1, "filtrate_recycle.fraction": 0.5},
    "water 5:1": {"pretreatment.water_ratio": 5},
}

# (지표, 모델 실행 조건, 문헌 하한, 상한, 출처)
LITERATURE = [
    ("phytate_conversion", {"enzyme.dose": 2500, "pretreatment.water_ratio": 5, "enzyme.pH": 5.5,
                            "enzyme.T": 50, "enzyme.time": 0.5},
     0.92, 0.92, "Liang 2009: 2.5 U/g, 1:5, pH 5.5, 50°C, 30 min → 피트산 92% 제거 (용해 포함)"),
    ("acidulation_P_recovery", {}, 0.94, 0.96, "EFMA BAT 4: 이수화물법 P2O5 회수 94-96%"),
    ("E1_steam_t_per_t_water", {}, 1.05, 1.33, "AIChE CEP 2018: 단효용 0.75-0.95 kg물/kg스팀"),
    ("E1_steam_t_per_t_P2O5", {}, 0.5, 2.2, "EFMA BAT 4: 농축 스팀 0.5-2.2 t/t P2O5"),
    ("E1_feed_P2O5_wt", {}, 0.26, 0.32, "EFMA BAT 4: 여과산 26-32% P2O5"),
]

# VAL-03 물성: (DB 상수, 문헌값, 출처)
PROPERTIES = [
    ("pKa1", 2.16, "CRC Handbook: H3PO4 pKa1"),
    ("pKa2", 7.21, "CRC Handbook: H3PO4 pKa2"),
    ("pKa3", 12.32, "CRC Handbook: H3PO4 pKa3"),
    ("HAP_logKsp", -58.33, "McDowell et al. 1977 J Res NBS 81A:273 권장값 4.7e-59 (25 °C)"),
]
PRECIP_PH = (6, 7, 8, 8.5, 9, 10)


def build_tables():
    base, sys, params = P.run()
    t = {"streams": P.stream_table(sys), "duties": P.duty_table(sys.units), "utilities": P.utility_table(sys, base)}

    rows = {name: P.run(o)[0] for name, o in SCENARIOS.items()}
    t["scenarios"] = pd.DataFrame(rows)
    feeds = {"rice_bran": rows["base"], **{fs: P.run(feedstock=fs)[0] for fs in P.FEEDSTOCKS}}
    t["feedstocks"] = pd.DataFrame(feeds)

    lit = []
    for key, overrides, lo, hi, ref in LITERATURE:
        m = P.run(overrides)[0][key] if overrides else base[key]
        lit.append({"지표": key, "모델": m, "문헌": f"{lo}" if lo == hi else f"{lo}-{hi}",
                    "판정": "범위 내" if lo * 0.95 <= m <= hi * 1.05 else "벗어남", "출처": ref})
    t["literature"] = pd.DataFrame(lit).set_index("지표")

    const = EQ.db_constants()
    t["properties"] = pd.DataFrame([{"상수": k, "DB": const[k], "문헌": ref_v, "차이": const[k] - ref_v, "출처": ref}
                                    for k, ref_v, ref in PROPERTIES]).set_index("상수")

    ext, v = sys.flowsheet.stream.extract, P.values(params)
    P_mM = ext.imol["H3PO4"] * 1e6 / ext.imass["Water"]           # mmol/kg 물
    Ca_mM = v["precipitation"]["Ca_P_ratio"] * P_mM
    t["precipitation"] = pd.DataFrame(
        {m: [EQ.removal(P_mM, Ca_mM, pH, v["precipitation"]["T"], m)["removal"] for pH in PRECIP_PH] for m in EQ.MINERALS},
        index=pd.Index(PRECIP_PH, name="pH"))
    precip_note = (f"추출액 P {P_mM:.1f} mM, Ca {Ca_mM:.1f} mM, {v['precipitation']['T']} °C. "
                   f"모델 가정 전환율 {v['precipitation']['conversion']}, 문헌 98.8-99% (pH 8.5, Ca/P 2; Warmadewanthi 2017)")

    t["parameters"] = pd.DataFrame([{"parameter": f"{g}.{k}", **{c: e.get(c) for c in ("value", "unit", "source", "ref", "range")}}
                                    for g, d in params.items() for k, e in d.items()]).set_index("parameter")
    summary = {
        "base": {k: (float(x) if not isinstance(x, str) else x) for k, x in base.items()},
        "parameter_count": {s: len(P.values_by_source(params, s)) for s in ("literature", "assumed", "measured")},
        "precipitation_note": precip_note,
    }
    return t, summary, sys


def show(t, summary):
    fmt = lambda df: df.map(lambda x: round(x, 3) if isinstance(x, float) else x).to_string()
    titles = {"streams": "스트림 (각 단계 투입·배출)", "duties": "단위공정 에너지", "utilities": "유틸리티·약품",
              "scenarios": "조건별 비교", "feedstocks": "원료별 비교", "literature": "문헌 비교 (VAL-02)",
              "properties": "물성 검증 (VAL-03): PHREEQC minteq.v4 vs 문헌", "precipitation": "침전 평형 P 제거율"}
    for k, title in titles.items():
        print(f"\n== {title}")
        if k == "precipitation":
            print("   " + summary["precipitation_note"])
        print(fmt(t[k]))
    params = t["parameters"]
    for source, label in (("literature", "문헌값"), ("assumed", "가정값")):
        sel = params[params["source"] == source]
        print(f"\n== {label} {len(sel)}개")
        for name, r in sel.iterrows():
            print(f"  {name} = {r['value']} {r['unit']}" + (f"  [{r['ref']}]" if r["ref"] else ""))


def save(t, summary, sys, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for k, df in t.items():
        df.to_csv(out / f"{k}.csv", encoding="utf-8-sig")  # Excel 한글 호환
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2))
    sys.diagram(file=str(out / "flowsheet"), format="png", display=False)
    print(f"\n저장: {out}/ ({', '.join(sorted(p.name for p in out.iterdir()))})")


if __name__ == "__main__":
    pd.set_option("display.width", 200)
    tables, summary, system = build_tables()
    system.diagram(file="flowsheet", format="png", display=False)
    print("공정도 저장: flowsheet.png")
    show(tables, summary)
    if "--out" in _sys.argv:
        save(tables, summary, system, _sys.argv[_sys.argv.index("--out") + 1])
