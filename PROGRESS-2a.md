# 2a 진행 상황 (끝나면 이 파일을 지우고 main에 합친다)

## 끝난 것
- 2a 브랜치 생성 (main 6fc0403 기준)
- relaxed scan 가능 여부 확인 → OCL 9.25.1(최신)로는 불가. 아래 "막힌 점" 참고.

## 다음에 할 것
- 사용자가 에너지 곡선 방식을 고르면 그에 맞춰 시작:
  1. Shape 탭 + 회전 가능한 C–C 결합 목록
  2. Newman 정렬 + 슬라이더 회전, staggered/eclipsed, anti/gauche
  3. 에너지 곡선 (SVG)
  4. 먼저 맞혀보기
  5. 각도 strain
  6. Shape 탭 떠나면 원래 모델로 복귀

## 막힌 점과 시도한 방법
- relaxed scan(이면각 고정 + 나머지 최적화)이 안 된다.
  - JS의 ForceFieldMMFF94는 minimise / getTotalEnergy / size만 노출. 생성자 options는 HashMap으로 무시된다.
  - Java의 AbstractForceField에 mFixedAtoms(원자 고정)가 있지만 JS 빌드에는 설정하는 코드가 없다(읽기만 남음).
    배포용 openchemlib.js는 내부 필드 이름이 난독화돼서 우회도 불가.
  - "짧게 minimise → 이면각 다시 맞추기" 반복도 불가: minimise는 수렴(반환 0)했을 때만 좌표를 분자에 쓴다.
    maxIts에 걸리면(반환 1) 좌표를 버린다.
  - RDKit.js(MinimalLib)에는 force field 함수가 없다.
  - 우리 최적화기 + OCL 에너지: 좌표마다 ForceFieldMMFF94를 새로 만들어야 해서(약 2–3 ms) 수치 미분이면 결합 하나에 몇 분. 불가.
- 고정 회전(rigid scan) 값 (MMFF94s+, kJ/mol, 각 분자의 최소 에너지 기준):
  - 에테인 장벽 13.5 (교과서 ~12) → 괜찮음
  - 뷰테인 gauche−anti 8.4 (교과서 ~3.8), CH3/CH3 가림형(0°) 44.3 (~19), H/CH3 가림형(120°) 17.3 (~16)
  - 자유 최적화한 극소점끼리는 맞음: 뷰테인 gauche −65.3°, gauche−anti 3.27
  → 고정 회전은 뷰테인에서 교과서 값과 크게 달라서 그대로 곡선을 낼 수 없다.

## 확인한 라이브러리 동작
- OCL minimise(): 반환 0 = 수렴, 1 = maxIts 도달, 2 = line search 실패. 0일 때만 분자 좌표가 바뀐다.
  build3D는 반환값을 보지 않지만, 테스트 분자(에테인~콜레스테롤, 사이클로알케인 3–6원)는 전부 0이었다.
- OCL MMFF 에너지 단위는 kcal/mol (×4.184 = kJ/mol).
- getCompactCopy()한 분자는 getAllConnAtoms 전에 ensureHelperArrays(cHelperNeighbours)가 필요하다.
- ForceFieldMMFF94 생성(타입 지정 포함)은 뷰테인에서 약 2–3 ms.
