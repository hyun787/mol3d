# 3b progress (delete before merging to main)

## Done
- 1. Isomer relationship ("Compare two structures", Stereo tab): `relation()` in chem.js, worker `{compare}` message, second structure as model 1 to the right (never with the mirror image), tests vs textbook answers and RDKit.

## Next
- 2. Meso mirror plane. 3. Fischer projection. 4. Guess first (relation, Fischer). Then README, merge, deploy, live check.

## Decisions
- Open stereo ("not set in the name") is checked only once the connectivity matches: constitutional isomers and non-isomers are told without stereo (1-butanol vs 2-butanol → constitutional isomers).

## Library behavior checked
- OCL `stripStereoInformation()` on a 3D copy does not remove stereo from `getIDCode()`: the code reads it back from the coordinates. `Molecule.fromIDCode(code)` then strip works (R = S = unset 2-butanol, E = Z).
- OCL parity from SMILES: 3 = left open, 1/2 = set, 0 = none. A C=C in a ring under 8 atoms gets parity 2 (never 3).
- `C[C@H]1C[C@@H]1C` is trans-1,2-dimethylcyclopropane; cis is `C[C@H]1C[C@H]1C`.
- OPSIN cannot be reached from the sandbox (and WebFetch needs the user's approval), so what names resolve to is checked on the live site.
