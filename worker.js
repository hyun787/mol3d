// Runs the slow chemistry (3D build, group search, ring flips) off the main thread so the page stays responsive.
import OCL from 'https://cdn.jsdelivr.net/npm/openchemlib@9.25.1/dist/openchemlib.js';
import { build3D, findGroups, atomInfo, kitParts, stereo, superimposable, newmanBonds, ringAngles, chairs, flipChair } from './chem.js';

const ready = OCL.Resources.registerFromUrl('https://cdn.jsdelivr.net/npm/openchemlib@9.25.1/dist/resources.json');
let last = null; // the last molecule built: a ring flip ({flip: chair index}) works on it

onmessage = async ({ data: { id, smiles, flip } }) => {
  try {
    await ready;
    if (flip !== undefined) { // only on request: the conformer search takes a second or more
      const f = flipChair(OCL, last.mol, last.chairs[flip].atoms);
      postMessage({ id, flip: f && { xyz: f.xyz, subs: f.subs, dE: f.dE } });
      return;
    }
    const { mol, stripped, minimized, given } = build3D(OCL, smiles);
    const atoms = atomInfo(mol);
    last = { mol, chairs: chairs(mol, atoms) };
    let charge = 0;
    for (let i = 0; i < mol.getAllAtoms(); i++) charge += mol.getAtomCharge(i);
    postMessage({
      id, stripped, minimized, charge,
      molfile: mol.toMolfile(),
      formula: mol.getMolecularFormula().formula,
      groups: findGroups(OCL, mol).map(({ group, matches, atoms }) => ({ id: group.id, count: matches.length, atoms })),
      atoms,
      kit: kitParts(mol),
      stereo: stereo(OCL, mol, given),
      superimposable: superimposable(mol),
      bonds: newmanBonds(OCL, mol, atoms),
      rings: ringAngles(mol, atoms),
      chairs: last.chairs,
    });
  } catch (e) {
    postMessage({ id, error: e?.message || String(e) });
  }
};
