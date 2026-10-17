# 3b progress (delete before merging to main)

## Done
- 1. Isomer relationship ("Compare two structures", Stereo tab): `relation()` in chem.js, worker `{compare}` message, second structure as model 1 to the right (never with the mirror image), tests vs textbook answers and RDKit.
- 2. Meso mirror plane: `isMeso()`, `mirrorPlane()` in chem.js (partners by the isotope substitution test; plane = bisector of two swapped stereocenters; 0.3 Å, one-to-one, no 1,4+ pair under 1.5 Å). Search order: model as built → turned about single bonds (inner chain bonds on a 5° grid, then coordinate descent, heavy atoms first) → generator seeds 1–20 (MMFF94s+, converged, same stereo), moved onto the model. Worker `{plane}`; button "Mirror plane" beside Mirror image; the model slides to the turned shape and the view turns to look along the plane.
- 3. Fischer projection: `fischer()` in chem.js (longest chain of non-ring carbons holding every chain stereocenter; more-oxidized end on top, then IUPAC-style locants; left/right from det(U − c, D − c, A − c) on the 3D model; pure alkyl tails written as one label; D/L for open-chain sugars by the lowest center and for α-amino acids by the α-carbon; meso drawings get a dashed mirror line). SVG in app.js with R/S tags; "Turn to the Fischer view" for one center (view only). Tests: textbook layouts, R/S read off the drawing = chem.js R/S, same drawing for another conformation.

## Next
- README, merge, deploy, live check (names through OPSIN/PubChem on the real site).

## Decisions
- Open stereo ("not set in the name") is checked only once the connectivity matches: constitutional isomers and non-isomers are told without stereo (1-butanol vs 2-butanol → constitutional isomers).
- The mirror plane is only drawn on coordinates that pass the check; cis-1,2-dimethylcyclohexane (chairs are chiral, ring bonds do not turn) gets the "another conformation" line plus a flat-ring note.
- D/L for amino acids uses the α-carbon, not the lowest stereocenter (the brief said lowest): L-threonine and L-isoleucine are L by C2. Sugars use the lowest stereocenter.

## Library behavior checked
- meso-tartaric acid built by MMFF94s+ is exactly centrosymmetric: one 180° turn of the central bond gives an exact plane (0.00 Å). meso-2,5-hexanediol needs the inner bond set first (coordinate descent alone stuck at 2.8 Å).
- OCL `stripStereoInformation()` on a 3D copy does not remove stereo from `getIDCode()`: the code reads it back from the coordinates. `Molecule.fromIDCode(code)` then strip works (R = S = unset 2-butanol, E = Z).
- OCL parity from SMILES: 3 = left open, 1/2 = set, 0 = none. A C=C in a ring under 8 atoms gets parity 2 (never 3).
- `C[C@H]1C[C@@H]1C` is trans-1,2-dimethylcyclopropane; cis is `C[C@H]1C[C@H]1C`.
- OPSIN cannot be reached from the sandbox (and WebFetch needs the user's approval), so what names resolve to is checked on the live site.
