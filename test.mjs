// Self-check: `npm install && node test.mjs`. Builds each molecule in 3D and checks the functional groups found.
import OCL from 'openchemlib';
import initRDKit from '@rdkit/rdkit';
import { build3D, findGroups, atomInfo, kitParts, stereo } from './chem.js';
import { resolveName } from './lookup.js';

OCL.Resources.registerFromNodejs();

// [name, SMILES, expected group ids with counts]
const CASES = [
  ['ethanol', 'CCO', 'alcohol'],
  ['acetic acid', 'CC(=O)O', 'acid'],
  ['aspirin', 'CC(=O)OC1=CC=CC=C1C(=O)O', 'aromatic ester acid'],
  ['benzaldehyde', 'O=Cc1ccccc1', 'aromatic aldehyde'],
  ['acetone', 'CC(C)=O', 'ketone'],
  ['ethyl acetate', 'CCOC(C)=O', 'ester'],
  ['acetamide', 'CC(N)=O', 'amide'],
  ['acetyl chloride', 'CC(Cl)=O', 'acylhalide'],
  ['acetic anhydride', 'CC(=O)OC(C)=O', 'anhydride'],
  ['aniline', 'Nc1ccccc1', 'aromatic amine'],
  ['triethylamine', 'CCN(CC)CC', 'amine'],
  ['phenol', 'Oc1ccccc1', 'aromatic phenol'],
  ['diethyl ether', 'CCOCC', 'ether'],
  ['ethanethiol', 'CCS', 'thiol'],
  ['dimethyl sulfide', 'CSC', 'sulfide'],
  ['nitrobenzene', 'O=[N+]([O-])c1ccccc1', 'aromatic nitro'],
  ['acetonitrile', 'CC#N', 'nitrile'],
  ['chlorobenzene', 'Clc1ccccc1', 'aromatic arylhalide'],
  ['2-chloropropane', 'CC(C)Cl', 'alkylhalide'],
  ['propene', 'C=CC', 'alkene'],
  ['propyne', 'C#CC', 'alkyne'],
  ['formaldehyde', 'C=O', 'aldehyde'],
  ['methyl formate', 'COC=O', 'ester'],
  ['formic acid', 'OC=O', 'acid'],
  ['N-methylethanimine', 'CC=NC', 'imine'],
  ['acetaldehyde oxime', 'CC=NO', ''],
  ['pyridine', 'c1ccncc1', 'aromatic'],
  ['pyrrole', 'c1cc[nH]c1', 'aromatic'],
  ['furan', 'c1ccoc1', 'aromatic'],
  ['thiophene', 'c1ccsc1', 'aromatic'],
  ['L-alanine', 'C[C@@H](C(=O)O)N', 'amine acid'],
  ['paracetamol', 'CC(=O)NC1=CC=C(C=C1)O', 'aromatic phenol amide'],
  ['ibuprofen', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', 'aromatic acid'],
  ['methyl acrylate', 'COC(=O)C=C', 'alkene ester'],
  ['acetophenone', 'CC(=O)c1ccccc1', 'aromatic ketone'],
  ['benzoyl chloride', 'ClC(=O)c1ccccc1', 'aromatic acylhalide'],
  ['methyl salicylate', 'COC(=O)C1=CC=CC=C1O', 'aromatic phenol ester'],
  ['dopamine', 'C1=CC(=C(C=C1CCN)O)O', 'aromatic phenol*2 amine'],
  ['oxirane', 'C1CO1', 'ether'],
  ['urea', 'NC(N)=O', 'amide*2'],
  ['sodium acetate', 'CC(=O)[O-].[Na+]', ''],
  ['2-aminopyridine', 'C1=CC=NC(=C1)N', 'aromatic amine'],
  ['imidazole', 'C1=CN=CN1', 'aromatic'],
  ['adenine', 'C1=NC2=NC=NC(=C2N1)N', 'aromatic*2 amine'],
  ['L-histidine', 'C1=C(NC=N1)C[C@@H](C(=O)O)N', 'aromatic amine acid'],
  ['cyclohexanone', 'O=C1CCCCC1', 'ketone'],
  ['indole', 'c1ccc2[nH]ccc2c1', 'aromatic*2'],
  ['L-tryptophan', 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N', 'aromatic*2 amine acid'],
  ['naphthalene', 'c1ccc2ccccc2c1', 'aromatic*2'],
  ['cholesterol', 'C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C', 'alkene alcohol'],
  ['caffeine', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'aromatic amide*3'],
  ['glucose', 'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O', 'alcohol*5 ether'],
];

let fail = 0;
for (const [name, smi, want] of CASES) {
  const t = performance.now();
  const { mol, minimized } = build3D(OCL, smi);
  const ms = (performance.now() - t).toFixed(0);
  const got = findGroups(OCL, mol).map(({ group, matches }) => group.id + (matches.length > 1 ? '*' + matches.length : '')).sort().join(' ');
  const exp = want.split(' ').filter(Boolean).sort().join(' ');
  const ok = got === exp && minimized;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', name.padEnd(20), `${ms}ms`.padStart(7), got, ok ? '' : `(want: ${exp}${minimized ? '' : ', not minimized'})`);
}
// Highlight set for ethanol's alcohol: O, its carbon, and the O–H hydrogen (not the C–H ones).
const eth = findGroups(OCL, build3D(OCL, 'CCO').mol)[0].atoms.length;
if (eth !== 3) { fail++; console.log('FAIL ethanol highlight atoms:', eth, '(want 3)'); }

// ---- Structure tab: textbook hybridization and lone pairs, heavy atoms in SMILES order ----
const HYB = [
  ['ethanol', 'CCO', 'sp3 sp3 sp3', '0 0 2'],
  ['acetic acid', 'CC(=O)O', 'sp3 sp2 sp2 sp2', '0 0 2 2'], // OH oxygen: lone pair in resonance with C=O
  ['acetamide', 'CC(N)=O', 'sp3 sp2 sp2 sp2', '0 0 1 2'],
  ['aniline', 'Nc1ccccc1', 'sp2 sp2 sp2 sp2 sp2 sp2 sp2', '1 0 0 0 0 0 0'],
  ['cis-2-butene', 'C/C=C\\C', 'sp3 sp2 sp2 sp3', '0 0 0 0'],
  ['cyclohexane', 'C1CCCCC1', 'sp3 sp3 sp3 sp3 sp3 sp3', '0 0 0 0 0 0'],
  ['benzene', 'c1ccccc1', 'sp2 sp2 sp2 sp2 sp2 sp2', '0 0 0 0 0 0'],
  ['acetonitrile', 'CC#N', 'sp3 sp sp', '0 0 1'],
  ['cyclopropane', 'C1CC1', 'sp3 sp3 sp3', '0 0 0'],
  ['triethylamine', 'CCN(CC)CC', 'sp3 sp3 sp3 sp3 sp3 sp3 sp3', '0 0 1 0 0 0 0'],
  ['diethyl ether', 'CCOCC', 'sp3 sp3 sp3 sp3 sp3', '0 0 2 0 0'],
  ['acetone', 'CC(C)=O', 'sp3 sp2 sp3 sp2', '0 0 0 2'],
  ['pyridine', 'c1ccncc1', 'sp2 sp2 sp2 sp2 sp2 sp2', '0 0 0 1 0 0'],
  ['pyrrole', 'c1cc[nH]c1', 'sp2 sp2 sp2 sp2 sp2', '0 0 0 1 0'],
  ['furan', 'c1ccoc1', 'sp2 sp2 sp2 sp2 sp2', '0 0 0 2 0'],
  ['methyl acetate', 'COC(C)=O', 'sp3 sp2 sp2 sp3 sp2', '0 2 0 0 2'],
  ['allyl amine', 'C=CCN', 'sp2 sp2 sp3 sp3', '0 0 0 1'], // N is next to CH2, not to the π bond
  ['propyne', 'C#CC', 'sp sp sp3', '0 0 0'],
  ['allene', 'C=C=C', 'sp2 sp sp2', '0 0 0'],
  ['tert-butyl cation', 'C[C+](C)C', 'sp3 sp2 sp3 sp3', '0 0 0 0'],
  ['water', 'O', 'sp3', '2'],
  ['ammonium', '[NH4+]', 'sp3', '0'],
  ['2-chloropropane', 'CC(C)Cl', 'sp3 sp3 sp3 -', '0 0 0 3'], // no hybridization shown for Cl
  ['dimethyl sulfoxide', 'CS(C)=O', 'sp3 - sp3 sp2', '0 1 0 2'],
];
for (const [name, smi, hyb, lp] of HYB) {
  const info = atomInfo(build3D(OCL, smi).mol).filter(x => x.el !== 'H');
  const got = info.map(x => x.hyb || '-').join(' '), gotLp = info.map(x => x.lp).join(' ');
  const ok = got === hyb && gotLp === lp;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'hyb', name.padEnd(20), got, '| lp', gotLp, ok ? '' : `(want: ${hyb} | ${lp})`);
}
// Bond angles measured on the 3D model: [SMILES, atom, test]
const ANGLES = [
  ['c1ccccc1', 0, a => a.every(x => Math.abs(x - 120) < 1), 'all 120°'],
  ['CC#N', 1, a => Math.abs(a[0] - 180) < 2, '180°'],
  ['C1CC1', 0, a => a.some(x => Math.abs(x - 60) < 1), 'C–C–C 60°'],
  ['CC', 0, a => a.every(x => Math.abs(x - 109.5) < 3), 'near 109.5°'],
  ['CC(N)=O', 2, a => Math.abs(a.reduce((s, x) => s + x) - 360) < 3, 'flat N (angles sum to 360°)'],
];
for (const [smi, at, test, want] of ANGLES) {
  const a = atomInfo(build3D(OCL, smi).mol)[at].angles.map(x => x[2]);
  const ok = test(a);
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'angle', smi.padEnd(18), a.map(x => x.toFixed(1)).join(' '), ok ? '' : `(want: ${want})`);
}

// ---- Kit tab: Darling Molecular Visions parts (nonzero counts) and how many items need checking ----
const KIT = [
  ['ethanol', 'CCO', 'black:4 red:2 h:6 lp:2 redOne:1', 0],
  ['acetic acid', 'CC(=O)O', 'black:2 red:2 half:1 end:1 h:4 lp:2 redOne:1', 0],
  ['acetamide', 'CC(N)=O', 'black:2 blue:2 half:1 end:1 h:5 lp:1', 0],
  ['aniline', 'Nc1ccccc1', 'blue:2 dbl:3 end:6 h:7 lp:1', 0],
  ['cis-2-butene', 'C/C=C\\C', 'black:4 dbl:1 end:2 h:8', 0],
  ['trans-2-butene', 'C/C=C/C', 'black:4 dbl:1 end:2 h:8', 0],
  ['cyclohexane', 'C1CCCCC1', 'black:12 h:12', 0],
  ['benzene', 'c1ccccc1', 'dbl:3 end:6 h:6', 0],
  ['acetonitrile', 'CC#N', 'black:2 h:3', 1], // C≡N is not in the kit
  ['cyclopropane', 'C1CC1', 'bend:3 h:6', 0],
  ['propyne', 'C#CC', 'black:2 tri:1 h:4', 0],
  ['pyridine', 'c1ccncc1', 'blue:1 dbl:3 end:5 h:5', 0],
  ['N-methylethanimine', 'CC=NC', 'black:4 blue:1 dbl:1 end:1 h:7', 0],
  ['2-chloropropane', 'CC(C)Cl', 'black:6 h:7', 1],
  ['cyclobutanone', 'O=C1CCC1', 'half:1 bend:3 h:6', 1],
];
for (const [name, smi, want, nCheck] of KIT) {
  const { n, check } = kitParts(build3D(OCL, smi).mol);
  const got = Object.entries(n).filter(([, v]) => v).map(([k, v]) => `${k}:${v}`).join(' ');
  const ok = got === want && check.length === nCheck;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'kit', name.padEnd(20), got, check.map(c => `[${c.what}]`).join(' '), ok ? '' : `(want: ${want}, ${nCheck} to check)`);
}

// ---- Stereo tab: E/Z cross-checked against RDKit's CIP labels on the same SMILES ----
const RDKit = await initRDKit();
const EZ = ['C/C=C\\C', 'C/C=C/C', 'OC(=O)/C=C\\C(=O)O', 'OC(=O)/C=C/C(=O)O', 'C/C(F)=C(/Cl)Br', 'C/C=N/O', 'C/C=N\\O',
  'CC(C)=CCC/C(C)=C/C=O', 'CC/C(C)=C\\C', 'CC/C(C(C)C)=C/C', 'OC/C(C)=C/C', 'C/C=C/C=C\\C', 'C1CCC/C=C\\CC1', 'OC(=O)/C=C/c1ccccc1',
  'CC/C(=C(\\c1ccccc1)/c1ccc(OCCN(C)C)cc1)/c1ccccc1', 'F/C=C(\\Cl)Br', 'C/C=C/C(=O)O', 'BrC/C=C/CCl', 'N#C/C=C\\C#N'];
for (const smi of EZ) {
  const { mol, given } = build3D(OCL, smi);
  const rd = RDKit.get_mol(smi);
  const z = JSON.parse(rd.get_json()).molecules[0].atoms.map(a => a.z ?? 6);
  const sameOrder = z.every((x, i) => mol.getAtomicNo(i) === x);
  const want = JSON.parse(rd.get_stereo_tags()).CIP_bonds.map(([p, q, l]) => `${Math.min(p, q)}-${Math.max(p, q)}${l[1]}`).sort().join(' ');
  rd.delete();
  const rows = stereo(OCL, mol, given).filter(r => r.kind === 'bond');
  const got = rows.map(r => `${Math.min(r.atoms[0], r.atoms[1])}-${Math.max(r.atoms[0], r.atoms[1])}${r.ez}`).sort().join(' ');
  const ctOk = rows.every(r => !r.ct || (r.ct === 'cis') === (r.ez === 'Z')); // 1,2-disubstituted: cis ⇔ Z
  const ok = sameOrder && got === want && ctOk && rows.every(r => r.given);
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'E/Z', smi.padEnd(22), got, rows.map(r => r.ct).filter(Boolean).join(' '), ok ? '' : `(RDKit: ${want}${sameOrder ? '' : ', atom order differs'})`);
}
// Left open in the SMILES: the model picks one, and the row says it was not given.
const open = stereo(OCL, build3D(OCL, 'CC=CC').mol, build3D(OCL, 'CC=CC').given);
if (open.length !== 1 || open[0].given) { fail++; console.log('FAIL 2-butene without stereo should be marked open'); }

// Ring cis/trans, from the 3D faces. Check: with two identical groups at 1,2 or 1,3, cis is the meso form,
// so RDKit's CIP labels at the two ring atoms differ (R,S); trans has the same label twice (R,R or S,S).
const RING = ['C[C@H]1CCCC[C@H]1C', 'C[C@H]1CCCC[C@@H]1C', 'C[C@@H]1CCC[C@H](C)C1', 'C[C@@H]1CCC[C@@H](C)C1',
  'O[C@H]1CCC[C@@H]1O', 'O[C@H]1CCC[C@H]1O', 'C[C@H]1C[C@@H]1C', 'C[C@H]1C[C@H]1C', 'Cl[C@H]1CCCC[C@@H]1Cl', 'Cl[C@H]1CCCC[C@H]1Cl'];
for (const smi of RING) {
  const { mol, given } = build3D(OCL, smi);
  const rd = RDKit.get_mol(smi);
  const cip = JSON.parse(rd.get_stereo_tags()).CIP_atoms.map(([, l]) => l.toUpperCase());
  rd.delete();
  const want = cip.length === 2 ? (cip[0] === cip[1] ? 'trans' : 'cis') : '?';
  const rows = stereo(OCL, mol, given).filter(r => r.kind === 'ring');
  const ok = rows.length === 1 && rows[0].ct === want && rows[0].given;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'ring', smi.padEnd(22), rows.map(r => r.ct).join(' '), ok ? '' : `(want: ${want})`);
}
const openRing = stereo(OCL, build3D(OCL, 'CC1CCC(C)CC1').mol, build3D(OCL, 'CC1CCC(C)CC1').given);
if (openRing.length !== 1 || openRing[0].given) { fail++; console.log('FAIL 1,4-dimethylcyclohexane without stereo should be marked open'); }

// Name lookup: OPSIN rejects, PubChem is busy (503), then unreachable, then answers.
const replies = [new Response('{}', { status: 404 }), new Response('', { status: 503 }), 'throw',
  new Response(JSON.stringify({ PropertyTable: { Properties: [{ CID: 2244, SMILES: 'CC(=O)OC1=CC=CC=C1C(=O)O', Title: 'Aspirin' }] } }))];
globalThis.fetch = async () => { const x = replies.shift(); if (x === 'throw') throw new TypeError('network'); return x; };
let waits = 0;
const hit = await resolveName('aspirin', async () => { waits++; });
if (hit?.cid !== 2244 || waits !== 2) { fail++; console.log('FAIL PubChem retry:', hit, 'waits', waits); }
globalThis.fetch = async () => new Response('', { status: 503 });
const busy = await resolveName('aspirin', async () => {}).then(() => 'no error', e => e.message);
if (!busy.startsWith('PubChem is busy')) { fail++; console.log('FAIL PubChem gives up:', busy); }

console.log(fail ? `${fail} failed` : 'all passed');
process.exit(fail ? 1 : 0);
