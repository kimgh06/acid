# acid — 폐농산물 유래 인산 공정 시뮬레이션 (1차)

미강(쌀겨)에서 인을 추출해 인산(H₃PO₄)으로 전환·농축하는 공정을 BioSTEAM으로 계산한다.
요구사항은 [todo.md](todo.md).

```
미강 + 물 + 피타아제
  → R1 효소 추출 → S1 고액분리 ─→ 추출 잔사(폐기)
  → R2 석회 침전 → S2 여과 ────→ 여액(폐기)
  → R3 황산 전환 → S3 석고 여과 ─→ 석고 케이크(폐기)
  → E1 3중효용 농축 → 제품 인산 + 응축수
```

> ⚠️ 현재 모든 매개변수는 **가정값**이다. 결과 수치는 모델 동작 확인용이며, 문헌·실측값으로 교체 전까지 설계 근거로 쓰지 않는다.

## 1. 설치

필요: [uv](https://docs.astral.sh/uv/), Graphviz(공정도용)

```bash
brew install uv graphviz     # Ubuntu: apt install graphviz + uv 공식 설치 스크립트
uv sync                      # Python 3.12 + 의존성(biosteam, phreeqc, pandas, pyyaml) 설치
```

## 2. 실행

```bash
uv run main.py           # 기준 공정 + 조건별 비교 결과 출력, flowsheet.png 저장
uv run test_balance.py   # 수지·조건 변화 검증 (모두 "ok" 가 떠야 정상)
```

`main.py` 출력 순서:

| 섹션 | 내용 |
|---|---|
| 스트림 | 각 단계 투입·배출 스트림의 유량(kg/h), 인(P kg/h), H₃PO₄(kg/h), 온도 |
| 단위공정 에너지 | 공정별 가열·냉각 부하(MJ/h), 전력(kW) |
| 조건별 비교 | 시나리오별 결과 지표(아래 4절) |
| 가정값 목록 | `source: assumed` 인 매개변수 전체 |

## 3. 조건 바꾸기

### 3.1 기준값 수정 — `params.yaml`

모든 항목은 `{value, unit, source}` 형식이다. 값만 바꾸고 다시 실행하면 된다.

```yaml
enzyme:
  T:  {value: 55, unit: degC, source: assumed}
```

| 그룹 | 주요 항목 |
|---|---|
| `feed` | 미강 투입량, 수분, 총인, 유기인 비율, 단백질·지질·탄수화물 (회분 = 나머지) |
| `pretreatment` | 물/미강 비 |
| `enzyme` | 온도, pH, 반응시간, 효소 투입량, 동역학 상수(k_ref, 최적 T·pH 등) |
| `separation1` | 잔사 케이크 함액률, 고형분 회수율, 미분해 피트산의 고상 잔류율 |
| `precipitation` | 석회 과잉률, 침전 전환율, 온도, 시간, 케이크 함액률 |
| `acidulation` | 황산 과잉률·농도, 전환율, 온도, 시간, 케이크 함액률 |
| `concentration` | 목표 H₃PO₄ 농도(wt) |

효소 전환율 모델: `X = 1 - exp(-k·t)`, `k = k_ref · (dose/dose_ref) · f_T · f_pH`
(f_T, f_pH 는 최적점 기준 가우시안 보정).

### 3.2 실측값 반영

실측·문헌값을 확보하면 값과 함께 `source` 를 바꾼다. 가정값 목록에서 자동으로 빠진다.

```yaml
  total_P: {value: 0.0213, unit: wt (dry), source: measured}   # 출처: OO연구소 2027-01 분석
```

### 3.3 시나리오 비교 — `main.py`

`SCENARIOS` 에 `"그룹.항목": 값` 형태로 추가한다. 기준값 대비 바뀐 항목만 적는다.

```python
SCENARIOS = {
    "base": {},
    "pH 4.5, 3h": {"enzyme.pH": 4.5, "enzyme.time": 3.0},
    "탈지미강": {"feed.lipid": 0.02},
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

## 4. 결과 지표

| 키 | 의미 |
|---|---|
| `phytate_conversion` | 피트산 → 인산염 전환율 |
| `P_extraction` | 원료 총인 중 추출액(S1 액상)으로 나온 비율 |
| `P_recovery` | **원료 총인 대비 제품 인 회수율** |
| `H3PO4_kg_h`, `product_kg_h`, `H3PO4_wt` | 제품 인산량, 제품 총량, 농도 |
| `purity_dry` | 무수 기준 순도 (H₃PO₄ / 물 제외 질량) |
| `P_loss_<스트림>_kg_h` | 공정별 인 손실 (잔사, 여액, 석고, 응축수) |
| `water_kg_h`, `lime_kg_h`, `H2SO4_kg_h`, `phytase_kg_h` | 물·약품 사용량 |
| `heating_MJ_h`, `cooling_MJ_h`, `electricity_kW` | 에너지 소비량 |
| `gypsum_kg_h`, `inositol_kg_h` | 부산물 |
| `mass_balance_error` | 전체 질량수지 상대오차 (≈ 0 이어야 함) |

`P_recovery × 원료 P + Σ P_loss = 원료 P` 가 항상 성립한다 (테스트로 검증).

## 5. 오류 메시지

| 오류 | 원인 · 조치 |
|---|---|
| `BalanceError: <공정> mass/P 수지 불일치` | 해당 단위공정에서 질량 또는 인이 보존되지 않음. 반응식·분리 모델 수정 필요 |
| `ValueError: 원료 조성 합 > 100 %` | `feed` 의 단백질+지질+탄수화물+인 성분 합이 건물 기준 100 % 초과 |
| `ValueError: 불순물이 많아 목표 농도 … 도달 불가` | 물을 다 빼도 목표 농도 미달. `target_wt` 를 낮추거나 불순물 원인 확인 |

## 6. 파일

| 파일 | 내용 |
|---|---|
| `params.yaml` | 원료·공정 조건 (가정/실측 구분) |
| `species.py` | 물질 정의. 단백질·지질·탄수화물·회분·효소는 대체물질 |
| `process.py` | 단위공정, 시스템 구성, 수지 검사, 결과·표 함수 |
| `main.py` | 기준 결과·시나리오 비교 출력, 공정도 저장 |
| `test_balance.py` | 질량·인 원소수지, 목표 농도, 조건 변화 방향 검증 |

## 7. 1차 한계

- 효소 동역학은 가우시안 T·pH 보정 1차 모델 — 문헌·실험값으로 교체 필요
- pH는 입력값. 산·염기 소비, 이온 평형, 침전 용해도는 미계산 (PHREEQC 연동 예정, `phreeqc` 설치됨)
- 미분해 피트산은 침전 공정에서 Ca-피테이트로 빠지지 않고 여액으로 손실 처리
- 석고 케이크 세척 없음 → 케이크 함액 중 인산 손실이 큼
- 반응열: 피타아제 가수분해 = 0 가정 (DB의 피트산 생성열 오류로 보정), 황산 희석열 미포함
- 분리기 전력, 장치 비용 미포함
