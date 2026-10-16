# Mol3D

Type a compound name and see the molecule in 3D, with buttons for the functional groups it contains. A static site for students: no login, works on phones, tablets and laptops.

## How it works

1. **Name → structure.** The browser asks [OPSIN](https://www.ebi.ac.uk/opsin/) (EMBL-EBI) to parse the name as a systematic (IUPAC) name. If OPSIN cannot parse it, the browser asks [PubChem](https://pubchem.ncbi.nlm.nih.gov/) to look it up as a common or trade name. Both services allow direct calls from the browser (CORS), so there is no server.
2. **Structure → 3D.** [OpenChemLib](https://github.com/cheminfo/openchemlib-js) builds a conformer and relaxes it with the MMFF94s+ force field, in a web worker so the page stays responsive. For salts and mixtures only the largest part is shown.
3. **Functional groups.** `chem.js` holds 21 groups (SMARTS + a short explanation). The groups found in the molecule become buttons that highlight their atoms.
4. **Viewer.** [3Dmol.js](https://3dmol.org).

## Files

- `index.html`: layout and styles
- `app.js`: search, viewer, buttons
- `worker.js`: 3D build and group search (OpenChemLib)
- `chem.js`: functional-group list and the chemistry helpers
- `lookup.js`: name → structure (OPSIN, then PubChem with retry)
- `test.mjs`: self-check over ~50 molecules (`npm install && npm test`)

## Notes

- A link like `…/mol3d/#caffeine` opens that molecule directly.
- PubChem allows about 5 requests per second per source. A whole class behind one school IP can hit that together, so `lookup.js` retries a busy reply up to 3 times after a random wait, which spreads the requests out. Names OPSIN can parse never reach PubChem.
- Molecules with large fused ring systems (steroids, taxol) take several seconds to build the first time.
