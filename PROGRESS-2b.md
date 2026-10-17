# 2b progress (delete before merging to main)

## Done
- [x] 1. Chair detection + list (Shape tab group "Chair (axial / equatorial)"). `chem.js`: `isChair` (6 ring torsions alternate in sign, each ≥ 30°), `chairs()` (6-rings of sp³ atoms that are chairs), `flipBlock` reasons. Tests: `chair …` rows.
- [x] 2. Axial/equatorial: `chairSubs` (angle to ring normal, |cos| > 0.65 = axial). "ax"/"eq" labels on the model, panel list, side-on view (`alignChair`). Newman bond, chair and ring-angle rows deselect each other. Tests: `ax/eq …` (textbook cis/trans ↔ a/e table against stereo()'s RDKit-checked cis/trans), `1,3-diaxial …` (nearest-atom cross-check).

## Next
- 3. Ring flip button (worker `{flip: i}` → `flipChair`), animation + caveat line, block reasons in the panel.
- 4. More stable chair (ΔE, |ΔE| < 1 kJ/mol = same), `judge` (noVerdict) reasons. Tests for the 10 validation cases + A-value table.
- 5. 1,3-diaxial rods/labels + panel line (`diaxial()` already in chem.js and tested).
- 6. Guess first for stability.
- Finish: README, delete this file, merge, deploy check, live check.

## Blocked / limits found
- OCL ConformerGenerator keeps only ring-fragment shapes close in energy: the flipped chair appears only for small groups (CH₃, OH, NH₂, F/Cl/Br/I, CN, dimethyls, 1,1-dimethyl). Never for tBu, Et, iPr, Ph, OMe, COOH, vinyl, menthol, β-glucose (40 seeds × all 4 strategies, 60° steps too). Asked the user whether a "reverse-pucker start + full MMFF minimisation" is acceptable; until answered, the spec method stays ("flip not available").
- MMFF94s+ Δ(ax−eq) kJ/mol vs A-values (Wikipedia table, kcal ×4.184): CH₃ 5.75 (7.1), OH 1.35 (3.6), NH₂ 2.81 (6.7) right side; Cl −1.45 (1.8), F −1.54 (0.6), Br −0.04 (1.6), CN −1.44 (0.7) wrong side; I +3.10 (1.8). Dimethyls: cis-1,2 / trans-1,3 / cis-1,4 = 0.00; trans-1,2 +7.54, cis-1,3 21.25, trans-1,4 11.71 favour diequatorial.
- Decision: verdict only when every non-H group is alkyl, OH or NH₂ (`noVerdict` in chem.js); halogens and C≡N say MMFF disagrees; anything else says not checked.

## Library behaviour confirmed
- `getOneConformerAsMolecule` ≈ 30 ms; `initializeConformers` ≈ 600 ms per seed (ring fragment step), 16 ms with `use60degreeSteps` on simple rings. Both keep atom order.
- `ForceFieldMMFF94` object can be reused: `ff.minimise()` then `ff.getTotalEnergy()` gives the minimised energy.
- build3D's model for methylcyclohexane and cis-1,3-dimethyl comes out axial (seed 42) — the model is not the lowest chair.
