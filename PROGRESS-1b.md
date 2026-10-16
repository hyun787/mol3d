# Phase 1b progress (delete this file when 1b merges into main)

Work happens on branch `1b`. main stays untouched until 1b is done.

## Steps
1. [x] R/S labels (Stereo tab rows, highlight center + neighbors, "not set" note)
2. [x] CIP priorities 1–4 (chem.js) + labels on the model + "put #4 in back" rotation
3. [ ] Guess-first toggle (off by default): R / S / Hint
4. [ ] Mirror image beside the model
5. [ ] README, delete this file, merge into main, check the live site

## Library behavior checked
- OCL 9.25.1 does not expose CIP ranks. `getFinalRanks`/`getSymmetryRank` are canonical ranks, not CIP.
- RDKit.js 2026.9.1 `get_json()` has `cipRanks` in the rdkitRepresentation extension, but these come from the
  legacy ranking, so they are not used as ground truth. `get_stereo_tags()` gives (R)/(S) and lowercase (r)/(s).
- OCL on pseudoasymmetric centers (cis/trans-1,4-dimethylcyclohexane): `isAtomParityPseudo` true and
  `getAtomCIPParity` 0, while RDKit gives r/s. They do not agree, so r/s is not shown.
- OCL also marks every inositol ring carbon as pseudo (cip 0) where RDKit gives R/S (they need CIP rule 4): skipped.
- OCL labels ribaric acid C3 (true pseudoasymmetric, RDKit "s") as plain S and not pseudo. Our own ranking ties
  there, which is why a center needs a full cipOrder to be shown.
- OCL marks the penicillin ring N (3 neighbors) as a stereocenter; RDKit does not. Only 4-neighbor centers are shown.
- OCL puts an atom parity on the central allene carbon (2 neighbors): skipped (4-neighbor centers only, so
  sulfoxide S is skipped too).
- OCL sulfoxide S and quaternary N R/S match RDKit. `[2H]` stays at its SMILES index in OCL (not moved to the end).

## Done so far
- chem.js `cipOrder` (rules 1a + 2 on the hierarchical digraph, duplicates for double bonds and ring closures)
  and `handedness` (R/S from the 3D model and the order). `stereo()` lists a center only when cipOrder separates
  all 4 groups AND its R/S equals OCL's. Checked against RDKit on 35 molecules (incl. cholesterol, taxol, morphine,
  strychnine) and 15 hand-checked priority orders in test.mjs.
- Stereo tab: "Stereocenters (R/S)" rows (name lists groups in atom order, not priority order, so guess mode
  can reuse it), tap = highlight center + neighbors + R/S label on the model.

- Step 2 UI: selecting a center puts 1–4 on its groups and R/S beside it (3Dmol labels), lists the priorities,
  and "Turn so 4 points away" animates the view (quaternion math in app.js; camera is on +z, so away = −z).
  Checked in the local browser: after the turn, group 4 has the same screen x,y as the center and lies behind it.

## Next
Step 3: guess-first toggle (off by default). On: hide R/S and priorities; selected center shows R / S / Hint.

## Local browser check
The sandbox cannot reach the CDNs, but Playwright + Chromium work: a harness routes jsDelivr to node_modules /
a downloaded 3dmol package and answers OPSIN with canned SMILES (scratchpad, not in the repo).

## Blocked / tried
(nothing yet)
