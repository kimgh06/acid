"""기준 공정 결과 + 조건별 비교 (SIM-09, SIM-10). 실행: uv run main.py"""
import pandas as pd

import process as P

SCENARIOS = {
    "base": {},
    "pH 4.0": {"enzyme.pH": 4.0},
    "pH 6.0": {"enzyme.pH": 6.0},
    "T 40C": {"enzyme.T": 40},
    "T 60C": {"enzyme.T": 60},
    "time 1h": {"enzyme.time": 1.0},
    "dose x2": {"enzyme.dose": 2000},
    "lime x1.2": {"precipitation.lime_excess": 1.2},
    "water 5:1": {"pretreatment.water_ratio": 5},
}

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
print(pd.DataFrame(rows).round(3).to_string())

assumed = P.assumed_values(params)
print(f"\n== 가정값 {len(assumed)}개 (실측 확보 시 params.yaml source: measured 로 교체)")
print("\n".join("  " + a for a in assumed))
