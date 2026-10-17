// Runs the slow chemistry (3D build, group search, ring flips) off the main thread so the page stays responsive.
import OCL from 'https://cdn.jsdelivr.net/npm/openchemlib@9.25.1/dist/openchemlib.js';
import { build3D, findGroups, atomInfo, kitParts, stereo, superimposable, relation, newmanBonds, ringAngles, chairs, flipChair, polarity, aromaticRings, equivalentH } from './chem.js';

const ready = OCL.Resources.registerFromUrl('https://cdn.jsdelivr.net/npm/openchemlib@9.25.1/dist/resources.json');
let last = null; // the last molecule built: a ring flip ({flip: chair index}) and a comparison ({compare: SMILES}) work on it

onmessage = async ({ data: { id, smiles, flip, compare } }) => {
  try {
    await ready;
    if (flip !== undefined) { // only on request: the conformer search takes a second or more
      const f = flipChair(OCL, last.mol, last.chairs[flip].atoms);
      // toMolfile() writes y and z negated (a 180° turn about x), so turn the flipped chair the same way for the viewer
      postMessage({ id, flip: f && { xyz: f.xyz.map(([x, y, z]) => [x, -y, -z]), subs: f.subs, dE: f.dE } });
      return;
    }
    if (compare !== undefined) { // a second structure, built in 3D for the viewer and compared with the last one
      const b = build3D(OCL, compare);
      postMessage({ id, compare: { molfile: b.mol.toMolfile(), stripped: b.stripped, rel: relation(OCL, last, b) } });
      return;
    }
    const { mol, stripped, minimized, given } = build3D(OCL, smiles);
    const atoms = atomInfo(mol);
    last = { mol, given, chairs: chairs(mol, atoms) };
    const pol = polarity(mol);
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
      pol: { ...pol, net: [pol.net[0], -pol.net[1], -pol.net[2]] }, // into the molfile's frame, as for flips
      arom: aromaticRings(mol),
      nmr: equivalentH(mol),
    });
  } catch (e) {
    postMessage({ id, error: e?.message || String(e) });
  }
};
