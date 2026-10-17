# Mol3D

Type a compound name and see the molecule in 3D, with buttons for the functional groups it contains. A static site for students: no login, works on phones, tablets and laptops.

## How it works

1. **Name → structure.** The browser asks [OPSIN](https://www.ebi.ac.uk/opsin/) (EMBL-EBI) to parse the name as a systematic (IUPAC) name. If OPSIN cannot parse it, the browser asks [PubChem](https://pubchem.ncbi.nlm.nih.gov/) to look it up as a common or trade name. Both services allow direct calls from the browser (CORS), so there is no server.
2. **Structure → 3D.** [OpenChemLib](https://github.com/cheminfo/openchemlib-js) builds a conformer and relaxes it with the MMFF94s+ force field, in a web worker so the page stays responsive. For salts and mixtures only the largest part is shown.
3. **Tabs, one per concept.**
   - **Structure:** buttons for the 21 functional groups in `chem.js` (SMARTS + a short explanation). Tap an atom to see its element, hybridization (textbook rule: a lone pair next to a π bond counts as sp²), lone pairs, and each bond angle measured on the model next to the ideal 109.5°/120°/180°.
   - **Kit:** the parts list for the Darling *Molecular Visions* kit, following Dr. Dan Berger's [guide](https://homepages.bluffton.edu/~bergerd/Models/visions.html) (Bluffton University). Aromatic rings count as one Kekulé structure. Parts the guide does not cover (halogens, S, P, C≡N, …) are listed to check against the kit.
   - **Stereo:** R/S for each stereocenter, E/Z for each stereo double bond (CIP rules), cis/trans for 1,2-disubstituted C=C, and cis/trans between substituted ring atoms. When the source structure leaves one open, the page says the 3D build picked it.
     - Tap a stereocenter to see its CIP priorities (1–4) on the model, and turn the model so group 4 points away. The priorities come from `chem.js` (rules 1 and 2 on the hierarchical digraph, with duplicate atoms); a center is shown only when they give the same R/S as OpenChemLib. Pseudoasymmetric centers (r/s) and the rare centers that need rules 3–5 are left out.
     - **Guess first** (off by default) hides R/S and the priorities: pick R or S to see the answer, or tap Hint to see the priorities and turn group 4 away first.
     - **Mirror image** shows the reflected molecule beside the model, R and S swapped, and says whether the two are superimposable (enantiomers, or one meso compound).
   - **Shape:** conformations and ring strain.
     - **Newman:** each C–C single bond outside a ring with two sp³ carbons (one per symmetry class). Tap one to look straight down it, front carbon toward you; a slider turns the back carbon, and the page says staggered or eclipsed, plus anti or gauche when each carbon carries exactly one non-H group.
     - **Energy curve:** schematic, order only (no numbers). It follows the textbook rule: eclipsed above staggered, and each eclipsed or gauche pair of carbon groups raises the energy. It is drawn only when every group on the two carbons is H, an sp³ carbon, or a C=C or benzene-ring carbon, since polar groups (C=O, O, N, halogens) can change the order. A C=C or benzene-ring group counts when eclipsed but not when gauche: MMFF94s+ puts its gauche form only 0.3–0.9 kJ/mol above anti (propylbenzene, 1-pentene), against 3.3 kJ/mol for butane, and `test.mjs` checks this. OpenChemLib cannot hold a dihedral fixed while it relaxes the rest (a relaxed scan), and a rigid scan overstates butane's barriers by about 2×, so the curve does not use force-field energies. **Guess first** (off by default) hides the curve until the student turns to the form they think is least stable.
     - **Steric (van der Waals) strain:** on bonds that get a curve, a see-through pink rod with a "vdW strain" label joins each pair of non-H groups sitting gauche or closer (Newman angle 65° or less; a C=C or benzene-ring group only when eclipsed), and the panel names them (gauche, eclipsed, or how far apart). Pairs with H are left out, following the textbook split: CH₃/CH₃ gauche is steric strain, H/CH₃ eclipsed is mostly torsional.
     - **Chair (axial / equatorial):** each six-membered ring of sp³ atoms that is a chair in the model (its six ring torsions alternate in sign, each at least 30°). Tap one to see it from the side, with "ax"/"eq" beside each group that is not H (axial = bond within about 50° of the ring's axis). Newman bonds, chairs and ring-angle rows deselect each other.
       - **1,3-diaxial strain:** each axial group (not H) gets a pink rod and a "1,3-diaxial" label to every axial atom two ring atoms away on the same face, and the panel names the pairs.
       - **Ring flip:** computed only when pressed, in the worker. OpenChemLib's conformer generator makes 40 conformers (one per seed), MMFF94s+ relaxes each fully (only converged runs count), and a chair with every group's axial/equatorial swapped is the flipped chair. The lowest one that keeps each group on its face is moved onto the model's ring; no coordinates are made up. The slide between the two chairs is a straight line, and the page says the real path goes through half-chair and twist-boat forms. The generator never builds the axial chair for groups larger than CH₃ or OH (tert-butyl, ethyl, isopropyl, phenyl, menthol, glucose), so those say "flip not available". Only-H rings, bridged rings and fused rings are not flipped, with the reason shown (cis-fused rings like cis-decalin can flip, but only with the other ring, which the page does not show).
       - **More stable chair:** the lowest MMFF94s+ energy of each chair, compared: "This chair is more stable", or "Both chairs are about the same" under 1 kJ/mol. It is given only for alkyl groups, OH and NH₂, which `test.mjs` checks against measured A-values (MMFF94s+ Δ axial − equatorial: CH₃ 5.8, NH₂ 2.8, OH 1.4 kJ/mol, same order as the A-values). MMFF94s+ puts F, Cl, Br and C≡N no lower equatorial (Cl −1.5 kJ/mol against an A-value of about +2), so halogens and nitriles get no verdict, and other groups get none until checked. All six dimethylcyclohexanes come out as in the textbook. **Guess first** (off by default) hides the answer until the student picks a chair.
     - **Ring angles:** for each ring, the angles at its atoms measured on the model, compared with the ideal (109.5° sp³, 120° sp²), with a short note by ring size.
4. **Viewer.** [3Dmol.js](https://3dmol.org).

## Files

- `index.html`: layout and styles
- `app.js`: search, viewer, buttons
- `worker.js`: 3D build, group search and ring flips (OpenChemLib)
- `chem.js`: functional-group list and the chemistry helpers
- `lookup.js`: name → structure (OPSIN, then PubChem with retry)
- `test.mjs`: self-check (`npm install && npm test`): functional groups over ~50 molecules, hybridization, lone pairs, bond angles, kit counts, Newman bond lists, their order and steric pairs (cross-checked against the 3D model and the textbook strain increments), ring angles, chairs (which rings count, axial/equatorial against the textbook cis/trans table, 1,3-diaxial pairs against distances, flipped chairs as true MMFF94s+ minima with every group swapped, which chair is more stable for the textbook cases, MMFF94s+ against A-values), hand-checked CIP priorities, and E/Z, ring cis/trans, R/S and mirror-image checks against [RDKit](https://www.rdkit.org) (a test-only dependency)

## Notes

- A link like `…/mol3d/#caffeine` opens that molecule directly.
- PubChem allows about 5 requests per second per source. A whole class behind one school IP can hit that together, so `lookup.js` retries a busy reply up to 3 times after a random wait, which spreads the requests out. Names OPSIN can parse never reach PubChem.
- Molecules with large fused ring systems (steroids, taxol) take several seconds to build the first time.
- The 3D build is one conformer, not always the lowest: methylcyclohexane and cis-1,3-dimethylcyclohexane come out with CH₃ axial. Ring flip shows the other chair and says which is more stable. A ring flip takes about a second.
