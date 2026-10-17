# 2a 추가 작업: vdW(steric) strain 표시 (끝나면 이 파일을 지우고 main에 합친다)

## 스펙 (사용자와 합의)
- 대상: 곡선이 나오는 결합만 (levels !== null, H와 sp³ C만). Cl/O/F 등은 gauche 효과 때문에 표시하지 않음.
- 판정은 거리 말고 각도: 결합 양쪽의 H 아닌 치환기 쌍(X/X)의 Newman 각도(슬롯 + 슬라이더 φ, 0–180°로 접음)가
  65° 이하(gauche ±5°부터 가림형까지)면 steric strain. 슬라이더 중간 각도(예: 30°)에서도 생겼다 사라짐.
- chem.js: stericPairs(bond, phi). newmanBonds의 슬롯 정보 재사용.
- 3D: 쌍마다 반투명 핑크 원통(--cut, 투명도 약 0.5)을 회전된 좌표로 그리고 가운데에 "vdW strain" 라벨.
  패널 영어 한 줄: "Steric (van der Waals) strain: CH₃ and CH₃ are gauche" / "... are eclipsed".
- 가림형에서 줄이 시선 방향으로 겹치는 건 처리 안 함 (사용자가 3D를 돌려서 봄).
- 근거(McMurry): CH₃/CH₃ gauche = steric, CH₃/CH₃ 가림 = torsional + steric, H/CH₃ 가림 = 대부분 torsional.

## 끝난 것
- 2a 브랜치 = main a2d01e4에서 시작
- 1단계: chem.js stericPairs + newmanBonds에 xx(X/X 쌍과 슬롯 차이) 추가, levels도 xx로 계산.
  test.mjs: 뷰테인 각도별 쌍 수, 에테인 0쌍, 순서 있는 모든 결합에서 정지점별 쌍 수를 3D 좌표 독립 계산과 대조,
  순서 없는 결합(DCE, 에탄올)은 쌍 없음. 기준을 125°로 바꾸면 7개 실패하는 것 확인.

## 다음에 할 것
2. app.js 원통·라벨·패널 문구, 브라우저 확인
3. README 한 줄 → 이 파일 삭제 → main에 합쳐 푸시 → 실제 사이트 확인

## 막힌 점과 시도한 방법
- (없음)

## 확인한 라이브러리 동작
- 1차 2a 기록은 main의 커밋 기록과 README 참고 (relaxed scan 불가, OCL minimise는 수렴 시에만 좌표 기록 등).
- 브라우저 확인: 스크래치패드 pw/harness.mjs (page.route로 OCL/3Dmol/OPSIN 응답). 스크린샷은 저장소 밖에.
