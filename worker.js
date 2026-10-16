// Runs the slow chemistry (3D build, group search) off the main thread so the page stays responsive.
import OCL from 'https://cdn.jsdelivr.net/npm/openchemlib@9.25.1/dist/openchemlib.js';
import { build3D, findGroups, atomInfo, kitParts, stereo } from './chem.js';

const ready = OCL.Resources.registerFromUrl('https://cdn.jsdelivr.net/npm/openchemlib@9.25.1/dist/resources.json');

onmessage = async ({ data: { id, smiles } }) => {
  try {
    await ready;
    const { mol, stripped, minimized, given } = build3D(OCL, smiles);
    let charge = 0;
    for (let i = 0; i < mol.getAllAtoms(); i++) charge += mol.getAtomCharge(i);
    postMessage({
      id, stripped, minimized, charge,
      molfile: mol.toMolfile(),
      formula: mol.getMolecularFormula().formula,
      groups: findGroups(OCL, mol).map(({ group, matches, atoms }) => ({ id: group.id, count: matches.length, atoms })),
      atoms: atomInfo(mol),
      kit: kitParts(mol),
      stereo: stereo(OCL, mol, given),
    });
  } catch (e) {
    postMessage({ id, error: e?.message || String(e) });
  }
};
