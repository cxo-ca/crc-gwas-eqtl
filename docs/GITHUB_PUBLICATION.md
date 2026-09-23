# GitHub 게시 기록

이 문서는 로컬 연구 폴더와 GitHub 저장소의 범위를 설명합니다.

## 저장소에 포함하는 연구 산출물

- `R/`, `scripts/`, `tests/`, `config.R`, `run.R`: 분석·검증 코드
- `docs/`: 설계, 자료 출처, 정렬, 공위치 분석, 감사 기록
- `reports/`: 최종 HTML 보고서, 소개 자료, 발표 자료
- `results/public/`, `results/final/`, `results/demo/`: 후보표, 최종 분류, 예제 결과
- `results/coloc/20260923_094957_real/`: 실제 분석 요약·QC·민감도·구간 그림. 단, 아래의 대용량 비압축 파일은 제외
- `data/processed/`: 분석용으로 정리한 자료와 체크섬
- `logs/`: 실행 환경과 패키지 버전 기록

## 커밋하지 않는 파일

|경로|이유|복구·대체 방법|
|---|---|---|
|`.Rlib/`|컴퓨터별 R 패키지 설치본|`scripts/setup.R` 실행|
|`data/raw/summary_statistics/`|다시 내려받을 수 있는 공개 원자료·구간 캐시 약 433 MB|`docs/DATA_ACQUISITION.md`, `docs/SOURCES.md`, `data/processed/download_manifest.tsv` 참고|
|`results/coloc/*/harmonized.tsv`|단일 파일 약 122 MB로 GitHub의 일반 Git 제한 초과|동일 정렬 입력의 압축본 `data/processed/harmonization/shared.aligned.tsv.gz` 사용|
|`data/raw/private/`, `data/input/*.tsv`|사용자 제공 자료가 놓일 수 있는 비공개 입력 위치|예제 입력과 `docs/DATA_DICTIONARY.md` 참고|
|`reports/introduction/.chrome-print/`|PDF 출력용 임시 브라우저 프로필|`scripts/build_intro_report.py`로 다시 생성|
|`tests/tmp*/`, `__pycache__/`|테스트·실행 중 생성되는 임시 파일|테스트와 스크립트 실행으로 재생성|

## 재현성 기록

원본과 분석 산출물의 SHA-256, 실행 설정, 경고, 세션 정보는 `data/processed/`, `results/coloc/20260923_094957_real/`, `results/final/`, `logs/`에 보존합니다. 최종 결과의 수치 검증 상태는 `results/final/final_summary.json`과 `results/final/final_validation.json`에서 확인할 수 있습니다.

이 제외 정책은 연구 결과를 생략하기 위한 것이 아니라, 제3자 공개 원자료의 중복 업로드와 GitHub 파일 크기 제한을 피하면서 분석 내용과 재현 경로를 모두 공개하기 위한 것입니다.
