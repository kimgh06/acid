# acid — 폐농산물 유래 인산 공정 시뮬레이션 (1차)

미강(쌀겨)에서 인을 추출해 인산(H₃PO₄)으로 전환·농축하는 공정을 BioSTEAM으로 계산한다.
요구사항은 [todo.md](todo.md).

```
미강 + 물(HX1 에서 여액 폐열로 예열) + 피타아제
  → R1 효소 추출 → S1 고액분리 ─→ 추출 잔사(폐기)
  → R2 석회 침전 → S2 여과 ────→ 여액 → HX1 → 폐기
  → R3 황산 전환 → S3 석고 여과 → W1 석고 세척 ─→ 석고 케이크(폐기)
  → M1 (여과산 + 세척액) → E1 단효용 진공 농축 → 제품 인산 + 응축수
```

> ⚠️ 매개변수 40개 중 24개는 **문헌값**, 16개는 **가정값**이다(`uv run main.py` 끝에 목록 출력).
> 실측 데이터로 검증되기 전까지 결과는 설계 근거가 아니라 경향 비교용이다.

## 1. 설치

필요: [uv](https://docs.astral.sh/uv/), Graphviz(공정도용)

```bash
brew install uv graphviz     # Ubuntu: apt install graphviz + uv 공식 설치 스크립트
uv sync                      # Python 3.12 + 의존성(biosteam, phreeqc, pandas, pyyaml) 설치
```

## 2. 실행

```bash
uv run main.py           # 기준 공정 + 조건별 비교 결과 출력, flowsheet.png 저장
uv run main.py --out reports/data/NNN   # 모든 표를 CSV·JSON 으로 저장 (보고서 부록)
uv run test_balance.py   # 수지·조건 변화 검증 (모두 "ok" 가 떠야 정상)
uv run sensitivity.py    # 가정값 ±20% 민감도 → 실측 우선순위
uv run uncertainty.py    # range 항목 몬테카를로 → 결과 P5/P50/P95, 영향 인자
uv run equilibrium.py    # PHREEQC 물성 상수·침전 평형
```

`main.py` 출력 순서:

| 섹션 | 내용 |
|---|---|
| 스트림 | 각 단계 투입·배출 스트림의 유량(kg/h), 인(P kg/h), H₃PO₄(kg/h), 온도 |
| 단위공정 에너지 | 공정별 가열·냉각 부하(MJ/h), 전력(kW) |
| 조건별 비교 | 시나리오별 결과 지표(아래 4절) |
| 문헌 비교 | 모델 결과 vs 문헌 범위, "범위 내/벗어남" 판정 (VAL-02, VAL-03) |
| 문헌값·가정값 목록 | `source: literature` (출처 포함), `source: assumed` 매개변수 |

## 3. 조건 바꾸기

### 3.1 기준값 수정 — `params.yaml`

모든 항목은 `{value, unit, source[, ref]}` 형식이다. `source` 는 `literature`(ref 필수) / `assumed` / `measured`. 값만 바꾸고 다시 실행하면 된다.

```yaml
enzyme:
  T:  {value: 55, unit: degC, source: assumed}
```

| 그룹 | 주요 항목 |
|---|---|
| `feed` | 미강 투입량, 수분, 총인, 유기인 비율, 단백질·지질·탄수화물 (회분 = 나머지) |
| `pretreatment` | 물/미강 비 |
| `enzyme` | 온도, pH, 반응시간, 효소 투입량, 동역학 상수(k_ref, 최적 T·pH, 실활 등) |
| `phosphatase` | 2차 효소(IP1 분해) 사용 여부·전환율 (기준 미사용) |
| `filtrate_recycle` | 여액 재순환 분율 (기준 0) |
| `heat_recovery` | 여액 → 공정수 예열 사용 여부, 최소 접근온도차 |
| `separation1` | 잔사 케이크 함액률, 고형분 회수율, 미분해 피트산의 고상 잔류율 |
| `precipitation` | Ca/P 몰비, 침전 전환율, 온도, 시간, 케이크 함액률 |
| `acidulation` | 황산 과잉률·농도, 전환율, 온도, 시간, 케이크 함액률, 세척 회수율·세척수 |
| `concentration` | 목표 H₃PO₄ 농도(wt), 증발기 압력 |

효소 모델 (피트산 → IP1 + 5 H₃PO₄, 피타아제는 IP1에서 멈춤):
- 전환율 `X = min(동역학, 효소용량)`
- 동역학 `1 - exp(-k·t_eff)`, `k = k_ref·(dose/dose_ref)·f_T·f_pH` (f_T, f_pH 는 최적점 기준 가우시안)
- 효소용량: 방출 Pi ≤ `dose × activity_ratio × f_T × f_pH × t_eff` (1 U = 1 µmol Pi/min)
- 열실활: `T > T_inact` 이면 `kd = kd_slope·(T - T_inact)`, `t_eff = (1 - e^(-kd·t))/kd`

결과의 `enzyme_limited_by` 가 어느 쪽이 전환율을 제한했는지 알려준다.

### 3.2 실측값 반영

실측·문헌값을 확보하면 값과 함께 `source` 를 바꾼다. 가정값 목록에서 자동으로 빠진다.

```yaml
  total_P: {value: 0.0213, unit: wt (dry), source: measured, ref: "OO연구소 2027-01 분석"}
```

### 3.3 시나리오 비교 — `main.py`

`SCENARIOS` 에 `"그룹.항목": 값` 형태로 추가한다. 기준값 대비 바뀐 항목만 적는다.

```python
SCENARIOS = {
    "base": {},
    "pH 4.5, 3h": {"enzyme.pH": 4.5, "enzyme.time": 3.0},
    "탈지미강": {"feed.lipid": 0.02, "feed.total_P": 0.022},
}
```

### 3.4 Python에서 직접 호출

```python
import process as P

result, sys, params = P.run({"enzyme.T": 45, "concentration.target_wt": 0.85})
result["P_recovery"]          # 지표 dict
P.stream_table(sys)           # 스트림 표 (pandas)
P.duty_table(sys.units)       # 에너지 표
P.assumed_values(params)      # 가정값 목록
sys.diagram(format="png")     # 공정도
```

`run()` 에 넘긴 값은 `source: scenario` 로 표시되어 가정값 목록에서 제외된다.

### 3.5 원료 바꾸기

`params.yaml` `feedstocks` 에 프리셋을 두고 `P.run(feedstock="defatted_rice_bran")` 으로 실행한다. 프리셋 값이 `feed` 를 덮어쓴다.
새 원료는 같은 6개 항목(수분, 총인, 유기인 비율, 단백질, 지질, 탄수화물)을 `source`·`ref` 와 함께 추가한다.

## 4. 결과 지표

| 키 | 의미 |
|---|---|
| `phytate_conversion` | 피트산 → IP1 전환율 |
| `phytate_P_released` | 피트산 P 중 무기인으로 방출된 비율 (피타아제만: 최대 5/6, 산성 포스파타아제 사용 시 그 이상) |
| `enzyme_limited_by` | 효소 전환율 제한 요인: `kinetics` / `capacity` |
| `P_extraction` | 순 추출률: 원료 총인 중 추출 잔사로 빠지지 않은 비율 (재순환 P 는 제외) |
| `P_recovery` | **원료 총인 대비 제품 인 회수율** (제품 중 모든 형태의 P) |
| `product_K_wt`, `product_Mg_wt`, `Mg_to_TCP_cake` | 제품 K·Mg 농도, 석회 침전 시 Mg 공침 분율 |
| `H3PO4_P_recovery` | 원료 총인 대비 H₃PO₄ 형태로 회수된 P (제품 중 IP1 등 유기 P 제외) |
| `H3PO4_kg_h`, `product_kg_h`, `H3PO4_wt` | 제품 인산량, 제품 총량, 농도 |
| `purity_dry` | 무수 기준 순도 (H₃PO₄ / 물 제외 질량) |
| `P_loss_<스트림>_kg_h` | 공정별 인 손실 (잔사, 여액, 석고, 응축수) |
| `acidulation_P_recovery` | 산 전환·세척 단계 P 회수율 (산 / Ca-인산 케이크) |
| `water_kg_h`, `lime_kg_h`, `H2SO4_kg_h`, `phytase_kg_h` | 물·약품 사용량 |
| `heating_MJ_h`, `cooling_MJ_h`, `electricity_kW` | 에너지 소비량 |
| `E1_steam_t_per_t_P2O5`, `E1_steam_t_per_t_water`, `E1_feed_P2O5_wt` | 농축 스팀 원단위, 증발기 입구 농도 |
| `gypsum_kg_h`, `IP1_kg_h`, `inositol_kg_h` | 부산물 |
| `mass_balance_error` | 전체 질량수지 상대오차 (≈ 0 이어야 함) |

`P_recovery × 원료 P + Σ P_loss = 원료 P` 가 항상 성립한다 (테스트로 검증).

## 5. 오류 메시지

| 오류 | 원인 · 조치 |
|---|---|
| `BalanceError: <공정> mass/P 수지 불일치` | 해당 단위공정에서 질량 또는 인이 보존되지 않음. 반응식·분리 모델 수정 필요 |
| `ValueError: 원료 조성 합 > 100 %` | `feed` 의 단백질+지질+탄수화물+인 성분 합이 건물 기준 100 % 초과 |
| `ValueError: 불순물이 많아 목표 농도 … 도달 불가` | 물을 다 빼도 목표 농도 미달. `target_wt` 를 낮추거나 불순물 원인 확인 |
| `ValueError: W1_gypsum_wash: 세척수 부족` | `acidulation.wash_water` 를 늘린다 |

## 6. 파일

| 파일 | 내용 |
|---|---|
| `params.yaml` | 원료·공정 조건 (가정/실측 구분) |
| `species.py` | 물질 정의. 단백질·지질·탄수화물·회분·효소·IP1은 대체물질 |
| `process.py` | 단위공정(Tank, Separator, Washer), 시스템 구성, 수지 검사, 결과·표 함수 |
| `main.py` | 기준 결과·시나리오 비교·문헌 비교 출력, 공정도 저장 |
| `sensitivity.py` | 가정값 ±20% 민감도 분석 |
| `uncertainty.py` | 몬테카를로 불확도 분석 (`range` 항목) |
| `equilibrium.py` | PHREEQC 침전 평형·물성 상수 검증 (VAL-03) |
| `reports/` | 단계별 보고서 (작성 규칙: `AGENTS.md`) |
| `AGENTS.md` | 에이전트 작업 규칙 (`CLAUDE.md` 가 import) |
| `test_balance.py` | 질량·인 원소수지, 목표 농도, 조건 변화 방향, 효소 종점, 세척 효과 검증 |

## 7. 결과·검증 현황

최신 결과와 해석은 단계 보고서에 있다 → [reports/README.md](reports/README.md) (단계별 P 회수율 추이), 최신 데이터는 `reports/data/<최신 번호>/`.
이 README 에는 수치를 적지 않는다 (단계마다 바뀌어 낡는다).

- 문헌 비교 판정: 모델값이 문헌 범위의 ±5 % 안이면 "범위 내" (`main.py` `LITERATURE`)
- 불확도: `uv run uncertainty.py` — 기준값과 함께 P5–P95 구간을 보고한다 (SR-008: 기준값이 분포 상단)
- 실측 우선순위: SR-003(1인자 민감도), SR-008(동시 변동) 9절 권고

## 8. 한계

- 효소: 가우시안 T·pH 보정 + 용량 상한 + 선형 실활. `k_ref`, `activity_ratio`, `kd_slope`, `T_sigma` 는 가정값
- pH는 입력값. 산·염기 소비는 미계산. 이온 평형은 R2 Mg 공침에만 PHREEQC 연동 (SR-012), P 침전 전환율은 문헌값
- 금속: K·Mg 만 추적 (Ca·Fe 등 미반영). 공침 Mg 의 추가 황산 소요 미반영
- 침전물은 Ca₃(PO₄)₂ 로 가정. 문헌상 pH 8.5–10 에서는 하이드록시아파타이트(Ca/P 1.67)가 생성될 수 있음
- 미분해 피트산·IP1은 침전 공정에서 빠지지 않고 여액으로 손실 처리
- 석고 세척은 회수율 하나로 요약, 세척액은 반응기 재순환 대신 농축기로 보냄
- 증발기: 끓는점오름 미반영 (압력으로 운전온도만 맞춤)
- 반응열: 피타아제 가수분해 = 0 가정 (DB의 피트산 생성열 오류로 보정), 황산 희석열 미포함
- 열회수: 여액 → 공정수 1기만. 전열계수는 BioSTEAM 기본값(0.5 kW/m²·K), 파울링 미반영
- 불확도: 매개변수 독립·균등분포 가정, 운전 조건은 표본 대상 아님
- 원료: 전지미강·탈지미강만 (밀기울·포도박은 문헌 조사 필요), 모든 원료에 같은 효소 동역학
- 분리기 전력, 장치 비용 미포함
