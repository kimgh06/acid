# ACID-SR-009 결과 자동 저장

| 문서번호 | Rev | 상태 | 작성 | 검토 | 일자 | BOD 기준 |
|---|---|---|---|---|---|---|
| ACID-SR-009 | A | 초안 | Claude | | 2026. 11. 16. | ACID-BOD Rev E |

커밋: 이 보고서와 같은 커밋 (`git log -1 -- reports/009-*.md`)

## 1. 요약

가. 결론
  1) `uv run main.py --out reports/data/NNN` 이 보고서에 쓰는 모든 표(11개 파일)를 CSV·JSON·PNG 로 저장한다. 보고서 수치를 콘솔에서 손으로 옮기던 과정을 없앴다.
  2) 유틸리티·약품 원단위 계산을 임시 스크립트에서 `process.utility_table` 로 옮기고, 전력 원단위 단위 오류(kWh/kg 을 t 당으로 표기)를 테스트로 고정했다.
  3) 공정 결과 변화 없음. 이후 단계 보고서는 `reports/data/NNN/` 을 부록으로 갖는다 (AGENTS.md 절차 1번).
나. 핵심 지표 (SR-007 전지미강 기준과 동일, 출처 `reports/data/009/summary.json`)
| 지표 | SR-008 기준값 | 이번 |
|---|---|---|
| P 회수율 (원료 총인 대비) | 0.670 | 0.670 |
| H3PO4 생산 (kg/h, 75 wt%) | 34.3 | 34.3 |
| 스팀 원단위 (t/t P2O5) | 14.6 | 14.6 |
| 가열 / 냉각 (MJ/h) | 776 / 410 | 776 / 410 |

## 2. 목적·범위

가. 질문: 보고서 수치를 실행 결과에서 직접, 재현 가능하게 가져올 수 있는가?
나. 다루지 않는 것: 불확도 결과 저장(`uncertainty.py` 는 콘솔 출력 유지), 보고서 본문 자동 생성.

## 3. 기준 및 가정

가. 매개변수·모델 변경 없음. 매개변수 현황: 문헌 24 / 가정 18 / 실측 0
나. CSV 인코딩 UTF-8 BOM (Excel 에서 한글 깨짐 방지).

## 4. 공정 설명

변경 없음.

## 5. 물질·열수지

`reports/data/009/streams.csv`, `duties.csv` — SR-005 와 동일 (질량·P 100.0 %, 에너지 잔차 1.5×10⁻⁸ kJ/h).

## 6. 유틸리티·약품 (`reports/data/009/utilities.csv`)

| 항목 | 사용량 | 원단위 (/t P2O5) |
|---|---|---|
| 스팀 | 362 kg/h | 14.6 t |
| 냉각수 | 5390 kg/h | 217 t |
| 전력 | 5.72 kW | 230 kWh |
| 공정수 + 세척수 | 8082 kg/h | 326 t |
| 석회 / 황산 / 효소 | 55.9 / 75.4 / 0.5 kg/h | 2.25 / 3.04 / 0.020 t |

## 7. 검증

가. `uv run test_balance.py` → 12 / 12 통과 (추가: `test_utility_table` — 전력 kWh/t, 약품 t/t 환산)
나. 저장 파일 11개 확인: duties, feedstocks, literature, parameters, precipitation, properties, scenarios, streams, utilities (CSV), summary.json, flowsheet.png
다. 위 6절 값이 SR-005 6절(임시 스크립트 계산)과 일치 → 이관 정확.

## 8. 목표 달성도

| # | 완료 조건 | 판정 | 근거 |
|---|---|---|---|
| 1 | 미강 투입량·성분 설정 | 달성 | |
| 2 | 전 공정 계산 | 달성 | |
| 3 | 공정별 인 이동·손실·회수량 | 달성 | `streams.csv` P 열 |
| 4 | 조건 변경 결과 차이 | 달성 | `scenarios.csv`, `feedstocks.csv` |
| 5 | 물질·에너지 수지 검증 | 달성 | |
| 6 | 인산 생산량·농도 | 달성 | |
| 7 | 문헌·실험 비교 | 부분 | `literature.csv`, `properties.csv`, 실측 없음 |
| 8 | 가정·실측 구분 | 달성 | `parameters.csv` source·ref·range 열 |

## 9. 결론 및 권고

| # | 권고 | 우선순위 | 근거 | 담당 |
|---|---|---|---|---|
| 1 | 실측값 입력 시 `parameters.csv` 를 실험 기록과 대조하는 검토 절차 | 중 | 가정·실측 구분의 감사 추적 | |
| 2 | `uncertainty.py` 에도 `--out` 추가 | 하 | 현재 콘솔만 | |

## 10. 참고문헌

없음 (도구 변경).

## 부록

가. 데이터: `reports/data/009/`

## 개정 이력

| Rev | 일자 | 내용 |
|---|---|---|
| A | 2026. 11. 16. | 최초 작성 |

끝.
