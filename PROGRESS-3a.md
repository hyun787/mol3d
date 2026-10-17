# 3a progress (delete before merging to main)

## Done
- chem.js: polarity(), aromaticRings() + huckel(), equivalentH(); test.mjs cases for all three (all pass).
- build3D fix: OCL's ConformerGenerator(42) built cis-1,2-dichloroethene (also F/Br analogues, some seeds) as trans.
  build3D now checks the 3D parities against the ones the SMILES set and tries seeds 1–20 when they differ.
  No existing test molecule needed a retry, so earlier results are unchanged.

- UI in Structure tab: Polarity, Aromatic rings (lobes; ring turned edge-on), Equivalent H (letters + 8 colors),
  guess-first for polar/nonpolar and signal count. Checked in Playwright harness (scratchpad h.mjs): shapes and labels
  go when leaving the tab, no page errors.

## Next
- Browser check with Playwright harness, README, delete this file, merge to main, check deploy and live site.

## Decisions / library behavior found
- Polar bond cutoff 0.4, not 0.5: with 0.5, C–N (0.49) and C≡N would not count, so pyridine, acetonitrile and
  trimethylamine would come out nonpolar. C–H (0.35) and S–H (0.38) stay out.
- Net cutoff 0.1 (sum of Δ·unit vectors): nonpolar cases sum < 0.02 (BF3 0.056: MMFF has no B types, raw conformer),
  polar ones ≥ 0.41.
- No polarity verdict for ions or molecules with atoms outside H B C N O F Cl Br (C–I, C–S: EN difference misleads).
- OCL 9.25.1 has only cHelperSymmetrySimple (63) and cHelperSymmetryStereoHeterotopicity (95); getSymmetryRank is
  undefined for H atoms in all modes, so it cannot count H environments. Used the D-substitution test instead:
  idcode of the D compound and of its mirror image; same pair = equivalent (enantiotopic merged, diastereotopic kept).
- OCL does not mark 2-pyridone aromatic; rings with an exocyclic C=O are skipped.
- npm test takes about 2 minutes (ring-flip searches); use a 10-minute timeout.
