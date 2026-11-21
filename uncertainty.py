"""불확도 분석: params.yaml 의 range 항목을 동시에 균등 표본추출 → 결과 분포와 주요 영향 인자.

실행: uv run uncertainty.py [N=300] [feedstock]
ponytail: 독립 균등분포 몬테카를로. 상관·분포형 정보가 생기면 LHS/정규분포로.
"""
import random
import sys

import pandas as pd

import process as P

METRICS = ["P_recovery", "H3PO4_kg_h", "heating_MJ_h", "E1_steam_t_per_t_P2O5", "purity_dry"]
SEED = 42


def ranged(params):
    return {f"{g}.{k}": v["range"] for g, d in params.items() for k, v in d.items() if "range" in v}


def sample(n, feedstock=None, seed=SEED):
    rng = random.Random(seed)
    ranges = ranged(P.load_params(feedstock=feedstock))
    rows, rejected = [], 0
    while len(rows) < n:
        x = {k: rng.uniform(lo, hi) for k, (lo, hi) in ranges.items()}
        try:
            r = P.run(x, feedstock=feedstock)[0]
        except ValueError:  # 조성 합 > 100 % 등 물리적으로 불가능한 조합
            rejected += 1
            continue
        rows.append({**x, **{m: r[m] for m in METRICS}})
    return pd.DataFrame(rows), list(ranges), rejected


def summarize(df, inputs):
    pct = df[METRICS].quantile([0.05, 0.5, 0.95]).T
    pct.columns = ["P5", "P50", "P95"]
    corr = df[inputs + METRICS].corr(method="spearman").loc[inputs, METRICS]
    return pct, corr


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    feedstock = sys.argv[2] if len(sys.argv) > 2 else None
    df, inputs, rejected = sample(n, feedstock)
    pct, corr = summarize(df, inputs)
    base = P.run(feedstock=feedstock)[0]
    pct.insert(0, "base", [base[m] for m in METRICS])
    pd.set_option("display.width", 200)
    print(f"표본 {n}개 (불가능 조합 {rejected}개 제외), 입력 {len(inputs)}개, seed {SEED}\n")
    print("== 결과 분포")
    print(pct.round(3).to_string())
    for m in ("P_recovery", "E1_steam_t_per_t_P2O5"):
        top = corr[m].abs().sort_values(ascending=False).head(5).index
        print(f"\n== {m} 영향 상위 5 (Spearman ρ)")
        print(corr.loc[top, m].round(3).to_string())


if __name__ == "__main__":
    main()
