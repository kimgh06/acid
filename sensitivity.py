"""가정값 민감도: 각 가정값을 ±20% 바꿔 결과 변화를 비교 → 실측 우선순위. 실행: uv run sensitivity.py"""
import pandas as pd

import process as P

STEP = 0.2
METRICS = ["P_recovery", "heating_MJ_h", "E1_steam_t_per_t_P2O5"]
SKIP = {"feed.bran_flow", "enzyme.dose_ref", "heat_recovery.enabled", "filtrate_recycle.fraction", "phosphatase.enabled"}  # 규모·기준값·on/off: ±20% 무의미


def perturb(value, sign):
    new = value * (1 + sign * STEP)
    return min(new, 0.999) if 0 < value < 1 else new  # 분율은 1 미만 유지


def main():
    base, _, params = P.run()
    rows = []
    for g, d in params.items():
        for k, v in d.items():
            key = f"{g}.{k}"
            if v["source"] != "assumed" or key in SKIP:
                continue
            lo, hi = (P.run({key: perturb(v["value"], s)})[0] for s in (-1, 1))
            row = {"parameter": key, "value": v["value"]}
            for m in METRICS:
                row[f"{m} -20%"] = lo[m] - base[m]
                row[f"{m} +20%"] = hi[m] - base[m]
            row["|ΔP_recovery| max"] = max(abs(row["P_recovery -20%"]), abs(row["P_recovery +20%"]))
            rows.append(row)
    df = pd.DataFrame(rows).sort_values("|ΔP_recovery| max", ascending=False).set_index("parameter")
    pd.set_option("display.width", 250)
    print(f"기준: P_recovery {base['P_recovery']:.3f}, heating {base['heating_MJ_h']:.0f} MJ/h, "
          f"steam {base['E1_steam_t_per_t_P2O5']:.2f} t/t P2O5\n")
    print("변화량(Δ) = 섭동 결과 - 기준. P_recovery 영향 큰 순.")
    print(df.round(3).to_string())


if __name__ == "__main__":
    main()
