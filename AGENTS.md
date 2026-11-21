# AGENTS.md

폐농산물(미강) 유래 인산 공정 시뮬레이션. BioSTEAM + Python, uv 환경. 요구사항은 `todo.md`, 사용법은 `README.md`.

## 명령

```bash
uv sync                  # 환경
uv run main.py           # 기준 결과 · 시나리오 비교 · 문헌 비교 · 가정값 목록
uv run main.py --out reports/data/NNN   # 단계 보고서 부록 데이터 저장 (CSV·JSON·공정도)
uv run test_balance.py   # 검증 — 모든 줄이 "ok" 여야 함
uv run sensitivity.py    # 가정값 ±20% 민감도
uv run uncertainty.py    # 불확도 (P5/P50/P95)
```

## 규칙

- 매개변수는 `params.yaml` 에만 둔다. 모든 항목에 `source` (`literature` 는 `ref` 필수 / `assumed` / `measured`). 코드에 숫자를 하드코딩하지 않는다 — 결과의 가정·실측 구분(todo.md 완료조건 8)이 깨진다.
- 단위공정을 바꾸면 `check_balance` (질량·P 원소수지)가 통과해야 한다. 수지를 맞추려고 허용오차를 키우지 않는다.
- 모듈 이름을 `chemicals.py` 로 짓지 않는다 — thermo 가 쓰는 `chemicals` 패키지를 가려 import 가 깨진다. (물질 정의는 `species.py`)
- 결과 수치를 보고할 때는 그 세션에서 실제로 실행한 출력만 쓴다.

## 단계 보고서 (매 단계 필수)

양식은 `reports/TEMPLATE.md` (설계 기준서 BOD + 단계 보고서 SR). 단계 = 모델 구조·매개변수·결과가 바뀌는 작업 묶음. 단계를 커밋할 때 함께:

1. `uv run main.py --out reports/data/NNN` 으로 부록 데이터를 저장하고, 보고서 수치는 이 파일에서 옮긴다.
2. `reports/NNN-<slug>.md` 를 SR 양식으로 새로 쓴다 (NNN 은 이전 번호 + 1).
3. BOD 항목(원료·제품 사양, 모델 기준, 단위공정 구성)이 바뀌었으면 `reports/BOD.md` 를 개정하고 Rev 를 올린다.
4. `reports/README.md` 표에 한 줄 추가한다.
5. 발행된 보고서는 고치지 않는다 — 결과가 바뀌면 새 보고서에 "이전 대비" 로 적는다.
