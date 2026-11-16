"""기준 공정 결과 + 조건별 비교(SIM-09, SIM-10) + 문헌 비교(VAL-02, VAL-03). 실행: uv run main.py"""
import pandas as pd

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

pd.set_option("display.width", 200)

base, sys, params = P.run()
sys.diagram(file="flowsheet", format="png", display=False)
print("공정도 저장: flowsheet.png\n")
print("== 스트림 (각 단계 투입·배출)")
print(P.stream_table(sys).round(3).to_string())
print("\n== 단위공정 에너지")
print(P.duty_table(sys.units).round(1).to_string())

print("\n== 조건별 비교")
rows = {name: P.run(o)[0] for name, o in SCENARIOS.items()}
print(pd.DataFrame(rows).map(lambda x: round(x, 3) if isinstance(x, float) else x).to_string())

print("\n== 문헌 비교 (VAL-02, VAL-03)")
lit = []
for key, overrides, lo, hi, ref in LITERATURE:
    m = P.run(overrides)[0][key] if overrides else base[key]
    lit.append({"지표": key, "모델": round(m, 3), "문헌": f"{lo}" if lo == hi else f"{lo}-{hi}",
                "판정": "범위 내" if lo * 0.95 <= m <= hi * 1.05 else "벗어남", "출처": ref})
print(pd.DataFrame(lit).to_string(index=False))

for source, label in (("literature", "문헌값"), ("assumed", "가정값")):
    items = P.values_by_source(params, source)
    print(f"\n== {label} {len(items)}개")
    print("\n".join("  " + a for a in items))
