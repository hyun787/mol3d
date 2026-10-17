// Self-check: `npm install && node test.mjs`. Builds each molecule in 3D and checks the functional groups found.
import OCL from 'openchemlib';
import initRDKit from '@rdkit/rdkit';
import { build3D, findGroups, atomInfo, kitParts, stereo, cipOrder, superimposable, relation, isMeso, mirrorPlane, fischer, newmanBonds, conformer, ringAngles, stericPairs, chairs, diaxial, flipChair, moreStable, polarity, aromaticRings, huckel, equivalentH, halides, sn2, sn2Rate, sn1, sn1Rate, e2, doubleBonds, addition } from './chem.js';
const dihedralOf = (m, a, b, c, d) => { const p = i => [m.getAtomX(i), m.getAtomY(i), m.getAtomZ(i)], s = (u, v) => u.map((x, k) => x - v[k]), cr = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
  const b2 = s(p(c), p(b)), n1 = cr(s(p(b), p(a)), b2), n2 = cr(b2, s(p(d), p(c))), x = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2], y = cr(n1, n2).reduce((t, v, k) => t + v * b2[k], 0) / Math.hypot(...b2);
  return Math.atan2(y, x) * 180 / Math.PI; };
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
  'CC/C(=C(\\c1ccccc1)/c1ccc(OCCN(C)C)cc1)/c1ccccc1', 'F/C=C(\\Cl)Br', 'C/C=C/C(=O)O', 'BrC/C=C/CCl', 'N#C/C=C\\C#N',
  'C1CCC=CC1', 'CC1=CC[C@@H](CC1)C(=C)C']; // 6-ring C=C: no E/Z row
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

// ---- R/S: every labeled center checked against RDKit's CIP labels (atoms by SMILES index), plus the textbook
// answer where the name fixes it. A center is only labeled when our priorities give the same R/S as OCL on the
// 3D model, so a passing row also checks the priorities. `skip` = centers RDKit labels that we leave out on purpose.
const RS = [
  ['L-alanine', 'C[C@@H](C(=O)O)N', 'S'],
  ['(R)-2-butanol', 'C[C@@H](O)CC', 'R'],
  ['(S)-2-butanol', 'C[C@H](O)CC', 'S'],
  ['D-glyceraldehyde', 'C([C@H](C=O)O)O', 'R'],
  ['L-serine', 'C([C@@H](C(=O)O)N)O', 'S'],
  ['(R)-2-chlorobutane', 'CC[C@@H](C)Cl', 'R'],
  ['(S)-2-chlorobutane', 'CC[C@H](C)Cl', 'S'],
  ['L-lactic acid', 'C[C@@H](C(=O)O)O', 'S'],
  ['(1R,2S,5R)-menthol', 'C[C@@H]1CC[C@H]([C@@H](C1)O)C(C)C', 'R S R'], // C5, C2, C1 in SMILES order
  ['cholesterol', 'C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C', 'R R S R S S R S'], // C20 C17 C14 C13 C9 C8 C10 C3
  ['cis-2-methylcyclohexanol', 'C[C@H]1CCCC[C@H]1O', ''],
  ['trans-2-methylcyclohexanol', 'C[C@@H]1CCCC[C@H]1O', ''],
  ['(2R,3R)-tartaric acid', '[C@@H]([C@H](C(=O)O)O)(C(=O)O)O', 'R R'],
  ['meso-tartaric acid', 'O[C@H](C(=O)O)[C@@H](O)C(=O)O', ''],
  ['3-methylpent-1-ene', 'C=C[C@H](C)CC', ''],
  ['3,4-dimethylpent-1-ene', 'C=C[C@H](C)C(C)C', ''],
  ['ethynyl vs tert-butyl', 'C#C[C@](C)(O)C(C)(C)C', ''],
  ['phenyl vs tert-butyl', 'c1ccccc1[C@H](O)C(C)(C)C', ''],
  ['(S)-1-deuterioethanol', '[2H][C@@H](C)O', 'S'],
  ['beta-D-glucose', 'C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O', 'R S S R R'],
  ['(S)-ibuprofen', 'C[C@@H](C1=CC=C(C=C1)CC(C)C)C(=O)O', 'S'],
  ['(1R,2S)-ephedrine', 'C[C@@H]([C@@H](C1=CC=CC=C1)O)NC', 'S R'], // C2 then C1
  ['L-threonine', 'C[C@H]([C@@H](C(=O)O)N)O', 'R S'], // C3 then C2
  ['L-isoleucine', 'CC[C@H](C)[C@@H](C(=O)O)N', 'S S'],
  ['limonene', 'CC1=CC[C@@H](CC1)C(=C)C', ''],
  ['camphor', 'C[C@@]12CC[C@@H](C1(C)C)CC2=O', ''],
  ['(S)-nicotine', 'CN1CCC[C@H]1C2=CN=CC=C2', 'S'],
  ['morphine', 'CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O', ''],
  ['penicillin G', 'CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C', ''], // the ring N is not listed
  ['taxol', 'CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C', ''],
  ['strychnine', 'C1CN2CC3=CCO[C@H]4CC(=O)N5[C@H]6[C@H]4[C@H]3C[C@H]2[C@@]61C7=CC=CC=C75', ''],
  ['ribaric acid (C3 is r/s)', 'O[C@H](C(=O)O)[C@H](O)[C@@H](O)C(=O)O', '', [5]],
  ['an inositol (needs CIP rule 4)', 'O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O', '', [1, 2, 4, 6, 8, 10]],
  ['cis-1,4-dimethylcyclohexane (r/s)', 'C[C@H]1CC[C@@H](C)CC1', ''],
  ['penta-2,3-diene (axial)', 'C/C=C=C/C', ''],
];
for (const [name, smi, named, skip = []] of RS) {
  const { mol, given } = build3D(OCL, smi);
  const rd = RDKit.get_mol(smi);
  const want = JSON.parse(rd.get_stereo_tags()).CIP_atoms.filter(([i, l]) => /[RS]/.test(l) && !skip.includes(i)).map(([i, l]) => i + l[1]).join(' ');
  rd.delete();
  const rows = stereo(OCL, mol, given).filter(r => r.kind === 'center');
  const got = rows.map(r => r.atom + r.rs).join(' ');
  const ok = got === want && rows.every(r => r.given) && (!named || rows.map(r => r.rs).join(' ') === named);
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'R/S', name.padEnd(34), got || '(none)', ok ? '' : `(RDKit: ${want || 'none'}${named ? `, name: ${named}` : ''})`);
}
// Left open in the SMILES: still labeled (the 3D build picked one), and marked as not given.
const openRS = stereo(OCL, build3D(OCL, 'CC(O)CC').mol, build3D(OCL, 'CC(O)CC').given).filter(r => r.kind === 'center');
if (openRS.length !== 1 || openRS[0].given) { fail++; console.log('FAIL 2-butanol without stereo should be one open center'); }

// CIP priorities checked by hand: [SMILES, center, groups from 1 to 4 (atom index, or H)]
const PRIO = [
  ['C[C@@H](C(=O)O)N', 1, [5, 2, 0, 'H'], 'alanine: NH2 > COOH > CH3 > H'],
  ['C([C@H](C=O)O)O', 1, [4, 2, 0, 'H'], 'glyceraldehyde: OH > CHO > CH2OH > H'],
  ['C([C@@H](C(=O)O)N)O', 1, [5, 2, 0, 'H'], 'serine: NH2 > COOH > CH2OH > H'],
  ['CC[C@@H](C)Cl', 2, [4, 1, 3, 'H'], '2-chlorobutane: Cl > ethyl > methyl > H'],
  ['C[C@@H](C(=O)O)O', 1, [5, 2, 0, 'H'], 'lactic acid: OH > COOH > CH3 > H'],
  ['C=C[C@H](C)CC', 2, [1, 4, 3, 'H'], 'vinyl > ethyl (duplicate C)'],
  ['C=C[C@H](C)C(C)C', 2, [1, 4, 3, 'H'], 'vinyl > isopropyl (decided one sphere later)'],
  ['C#C[C@](C)(O)C(C)(C)C', 2, [4, 1, 5, 3], 'OH > ethynyl > tert-butyl > CH3'],
  ['c1ccccc1[C@H](O)C(C)(C)C', 6, [7, 5, 8, 'H'], 'OH > phenyl > tert-butyl > H'],
  ['O=C[C@H](Cl)C(O)C', 2, [3, 1, 4, 'H'], 'Cl > CHO {O,O,H} > CH(OH)CH3 {O,C,H} > H'],
  ['O[C@@H](CCBr)C(C)C', 1, [0, 5, 2, 'H'], 'isopropyl {C,C,H} > CH2CH2Br {C,H,H}: first point of difference'],
  ['[C@@H](O)(CBr)C(C)C', 0, [1, 2, 4, 'H'], 'CH2Br {Br,H,H} > isopropyl {C,C,H}'],
  ['C[C@H]1CCCC=C1', 1, [6, 2, 0, 'H'], '3-methylcyclohexene: ring CH= side > ring CH2 side'],
  ['C[C@H]1CCCCC1=O', 1, [6, 2, 0, 'H'], '2-methylcyclohexanone: C=O side > CH2 side'],
  ['[2H][C@@H](C)O', 1, [3, 2, 0, 'H'], 'OH > CH3 > D > H (rule 2: mass)'],
];
// First-sphere sets shown with the priorities: [SMILES, center, sets from priority 1 to 4]
const SETS = [
  ['C[C@@H](C(=O)O)N', 1, 'H,H | O,O,O | H,H,H | '], // alanine: N, COOH (C=O counts twice), CH3, H
  ['C([C@H](C=O)O)O', 1, 'H | O,O,H | O,H,H | '], // glyceraldehyde: OH, CHO, CH2OH, H
  ['C=C[C@H](C)CC', 2, 'C,C,H | C,H,H | H,H,H | '], // vinyl (duplicate C), ethyl, methyl, H
];
for (const [smi, c, want] of SETS) {
  const { mol, given } = build3D(OCL, smi);
  const got = stereo(OCL, mol, given).find(r => r.atom === c)?.sets.map(x => x.join()).join(' | ');
  if (got !== want) fail++;
  console.log(got === want ? 'ok  ' : 'FAIL', 'sets', smi.padEnd(22), got, got === want ? '' : `(want ${want})`);
}
for (const [smi, c, want, what] of PRIO) {
  const { mol } = build3D(OCL, smi);
  const got = (cipOrder(mol, c) || []).map(a => (a >= mol.getAtoms() ? 'H' : a));
  const ok = got.join() === want.join();
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'priority', what, ok ? '' : `(got ${got.join(' ')}, want ${want.join(' ')})`);
}

// Mirror image: superimposable or not. RDKit check: invert every @/@@ and compare canonical SMILES.
const MIRROR = [
  ['L-alanine', 'C[C@@H](C(=O)O)N'], ['(2R,3R)-tartaric acid', '[C@@H]([C@H](C(=O)O)O)(C(=O)O)O'],
  ['meso-tartaric acid', 'O[C@H](C(=O)O)[C@@H](O)C(=O)O'], ['meso-2,3-dibromobutane', 'C[C@@H](Br)[C@H](C)Br'],
  ['cis-1,2-dimethylcyclohexane', 'C[C@H]1CCCC[C@H]1C'], ['trans-1,2-dimethylcyclohexane', 'C[C@H]1CCCC[C@@H]1C'],
  ['ribaric acid', 'O[C@H](C(=O)O)[C@H](O)[C@@H](O)C(=O)O'], ['cholesterol', CASES.find(c => c[0] === 'cholesterol')[1]],
  ['(1R,2S)-ephedrine', 'C[C@@H]([C@@H](C1=CC=CC=C1)O)NC'], ['ethanol', 'CCO'], ['cis-2-butene', 'C/C=C\\C'],
];
for (const [name, smi] of MIRROR) {
  const can = s => { const m = RDKit.get_mol(s), c = m.get_smiles(); m.delete(); return c; };
  const want = can(smi) === can(smi.replace(/@@?/g, x => (x === '@' ? '@@' : '@')));
  const got = superimposable(build3D(OCL, smi).mol);
  if (got !== want) fail++;
  console.log(got === want ? 'ok  ' : 'FAIL', 'mirror', name.padEnd(30), got ? 'superimposable' : 'not superimposable', got === want ? '' : '(RDKit disagrees)');
}

// ---- How two structures are related: textbook answer, and RDKit for the stereo verdicts (canonical isomeric SMILES;
// the mirror image inverts every @/@@ and keeps / and \). Sugars in the open-chain form, checked by name via RDKit CIP below.
const SMI = {
  '(R)-2-butanol': 'C[C@@H](O)CC', '(S)-2-butanol': 'C[C@H](O)CC', '2-butanol': 'CC(O)CC', '1-butanol': 'CCCCO',
  '(2R,3R)-tartaric acid': '[C@@H]([C@H](C(=O)O)O)(C(=O)O)O', '(2S,3S)-tartaric acid': '[C@H]([C@@H](C(=O)O)O)(C(=O)O)O',
  'meso-tartaric acid': 'O[C@H](C(=O)O)[C@@H](O)C(=O)O', 'meso-tartaric acid (written the other way)': 'O[C@@H](C(=O)O)[C@H](O)C(=O)O',
  'cis-2-butene': 'C/C=C\\C', 'trans-2-butene': 'C/C=C/C', '2-butene': 'CC=CC',
  'D-glucose': 'O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO', 'D-galactose': 'O=C[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO',
  'L-glucose': 'O=C[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)CO',
  'beta-D-glucopyranose': 'C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O', 'beta-D-galactopyranose': 'C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O',
  'cis-1,2-dimethylcyclohexane': 'C[C@H]1CCCC[C@H]1C', 'cis-1,2-dimethylcyclohexane (mirror-written)': 'C[C@@H]1CCCC[C@@H]1C',
  'trans-1,2-dimethylcyclohexane': 'C[C@H]1CCCC[C@@H]1C',
  butane: 'CCCC', isobutane: 'CC(C)C', '1-propanol': 'CCCO', '2-propanol': 'CC(C)O', ethanol: 'CCO', 'dimethyl ether': 'COC', methanol: 'CO',
};
// The sugar SMILES carry the names' configurations (RDKit CIP, atoms by SMILES index).
for (const [name, want] of [['D-glucose', '2R 4S 6R 8R'], ['D-galactose', '2R 4S 6S 8R'], ['L-glucose', '2S 4R 6S 8S']]) {
  const rd = RDKit.get_mol(SMI[name]), got = JSON.parse(rd.get_stereo_tags()).CIP_atoms.map(([i, l]) => i + l[1]).join(' ');
  rd.delete();
  if (got !== want) { fail++; console.log('FAIL', name, 'SMILES has', got, `(want ${want})`); }
}
const RELATION = [
  ['(R)-2-butanol', '(S)-2-butanol', 'enantiomers'], ['(2R,3R)-tartaric acid', 'meso-tartaric acid', 'diastereomers'],
  ['meso-tartaric acid', 'meso-tartaric acid (written the other way)', 'same'], ['(2R,3R)-tartaric acid', '(2S,3S)-tartaric acid', 'enantiomers'],
  ['cis-2-butene', 'trans-2-butene', 'diastereomers'], ['D-glucose', 'D-galactose', 'diastereomers'], ['D-glucose', 'L-glucose', 'enantiomers'],
  ['beta-D-glucopyranose', 'beta-D-galactopyranose', 'diastereomers'],
  ['cis-1,2-dimethylcyclohexane', 'trans-1,2-dimethylcyclohexane', 'diastereomers'],
  ['cis-1,2-dimethylcyclohexane', 'cis-1,2-dimethylcyclohexane (mirror-written)', 'same'], // meso
  ['butane', 'isobutane', 'constitutional'], ['1-propanol', '2-propanol', 'constitutional'], ['ethanol', 'dimethyl ether', 'constitutional'],
  ['ethanol', 'methanol', 'not isomers'], ['(R)-2-butanol', '(R)-2-butanol', 'same'],
  ['2-butanol', '(R)-2-butanol', 'open'], ['cis-2-butene', '2-butene', 'open'], ['1-butanol', '2-butanol', 'constitutional'],
];
const rdCan = s => { const m = RDKit.get_mol(s), c = m.get_smiles(); m.delete(); return c; };
const rdMass = s => { const m = RDKit.get_mol(s), d = JSON.parse(m.get_descriptors()); m.delete(); return d.exactmw.toFixed(4); };
const rdRelation = (a, b) => {
  if (rdMass(a) !== rdMass(b)) return 'not isomers';
  const flat = s => rdCan(s.replace(/[@/\\]/g, ''));
  if (flat(a) !== flat(b)) return 'constitutional';
  if ([a, b].some(s => s === SMI['2-butanol'] || s === SMI['2-butene'])) return 'open'; // the structures with open stereo, listed by hand
  return rdCan(a) === rdCan(b) ? 'same' : rdCan(a) === rdCan(b.replace(/@@?/g, x => (x === '@' ? '@@' : '@'))) ? 'enantiomers' : 'diastereomers';
};
for (const [x, y, want] of RELATION) {
  const got = relation(OCL, build3D(OCL, SMI[x]), build3D(OCL, SMI[y]));
  const rd = rdRelation(SMI[x], SMI[y]), ok = got.kind === want && rd === want;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'relation', `${x} / ${y}`.padEnd(58), got.kind, ok ? '' : `(want ${want}, RDKit ${rd})`);
}
{ // two conformations of one compound: (R)-2-butanol with its C2–C3 bond turned 120°, a rigid turn as on the Shape tab
  const b = build3D(OCL, SMI['(R)-2-butanol']), other = b.mol.getCompactCopy(), t0 = dihedralOf(b.mol, 0, 1, 3, 4);
  other.ensureHelperArrays(OCL.Molecule.cHelperNeighbours);
  const P = i => [other.getAtomX(i), other.getAtomY(i), other.getAtomZ(i)], c = P(3), k = P(3).map((x, j) => x - P(1)[j]), n = Math.hypot(...k);
  const side = new Set([3]); // atoms beyond C3, away from C2
  for (const a of side) for (let j = 0; j < other.getAllConnAtoms(a); j++) { const x = other.getConnAtom(a, j); if (x !== 1) side.add(x); }
  side.delete(3);
  const th = 2 * Math.PI / 3, u = k.map(x => x / n), cs = Math.cos(th), sn = Math.sin(th);
  for (const a of side) {
    const v = P(a).map((x, j) => x - c[j]), w = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]], d = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    const r = v.map((x, j) => x * cs + w[j] * sn + u[j] * d * (1 - cs) + c[j]);
    other.setAtomX(a, r[0]); other.setAtomY(a, r[1]); other.setAtomZ(a, r[2]);
  }
  const t = dihedralOf(other, 0, 1, 3, 4), got = relation(OCL, b, { mol: other, given: b.given }).kind;
  if (got !== 'same') { fail++; console.log('FAIL relation: two conformations of (R)-2-butanol should be the same compound, got', got); }
  else console.log('ok   relation: (R)-2-butanol in two conformations', `C–C–C–C ${t0.toFixed(0)}° and ${t.toFixed(0)}°`, 'same');
}

// ---- Meso mirror plane: [name, SMILES, 'plane' | 'none' (meso, but no plane on reachable shapes) | 'not meso']
// A found plane is checked on its own: reflecting every atom through it gives the same set of positions (each atom onto
// a different atom of its element, within 0.3 Å); bonds keep their model lengths; the stereo of the shape is the
// model's; and the plane has a stereocenter on each side.
const MESO = [
  ['meso-tartaric acid', SMI['meso-tartaric acid'], 'plane'], ['meso-2,3-dibromobutane', 'C[C@@H](Br)[C@H](C)Br', 'plane'],
  ['meso-2,4-pentanediol', 'C[C@@H](O)C[C@H](C)O', 'plane'], ['cis-1,2-dimethylcyclopropane', 'C[C@H]1C[C@H]1C', 'plane'],
  ['cis-1,3-dimethylcyclohexane', 'C[C@@H]1CCC[C@H](C)C1', 'plane'], ['ribaric acid (plane through the r/s center)', 'O[C@H](C(=O)O)[C@H](O)[C@@H](O)C(=O)O', 'plane'],
  ['meso-2,5-hexanediol', 'C[C@@H](O)CC[C@H](C)O', 'plane'],
  ['cis-1,2-dimethylcyclohexane (chairs are chiral)', SMI['cis-1,2-dimethylcyclohexane'], 'none'],
  ['(2R,3R)-tartaric acid', SMI['(2R,3R)-tartaric acid'], 'not meso'], ['(R)-2-butanol', SMI['(R)-2-butanol'], 'not meso'],
  ['cis-1,4-dimethylcyclohexane (no stereocenter)', 'C[C@H]1CC[C@@H](C)CC1', 'not meso'],
];
for (const [name, smi, want] of MESO) {
  const { mol } = build3D(OCL, smi), t = performance.now(), meso = isMeso(OCL, mol), r = mirrorPlane(OCL, mol), ms = (performance.now() - t).toFixed(0);
  const got = !meso ? 'not meso' : r ? 'plane' : 'none';
  let ok = got === want, why = '';
  if (r) {
    const n = mol.getAllAtoms(), el = a => mol.getAtomicNo(a), P = r.xyz, used = new Set();
    for (let a = 0; a < n; a++) {
      const d = 2 * ((P[a][0] - r.c[0]) * r.n[0] + (P[a][1] - r.c[1]) * r.n[1] + (P[a][2] - r.c[2]) * r.n[2]), q = P[a].map((x, k) => x - d * r.n[k]);
      const b = [...Array(n).keys()].filter(j => el(j) === el(a)).sort((i, j) => Math.hypot(...P[i].map((x, k) => x - q[k])) - Math.hypot(...P[j].map((x, k) => x - q[k])))[0];
      if (Math.hypot(...P[b].map((x, k) => x - q[k])) > 0.3) why ||= `atom ${a} reflects onto nothing`;
      used.add(b);
    }
    if (used.size !== n) why ||= 'two atoms reflect onto one';
    const len = (X, b) => Math.hypot(...X[mol.getBondAtom(0, b)].map((x, k) => x - X[mol.getBondAtom(1, b)][k]));
    const home = [...Array(n).keys()].map(a => [mol.getAtomX(a), mol.getAtomY(a), mol.getAtomZ(a)]);
    for (let b = 0; b < mol.getAllBonds(); b++) if (Math.abs(len(P, b) - len(home, b)) > (r.how === 'conformer' ? 0.05 : 0.001)) why ||= `bond ${b} length changed`;
    const m = mol.getCompactCopy(); P.forEach(([x, y, z], a) => { m.setAtomX(a, x); m.setAtomY(a, y); m.setAtomZ(a, z); });
    if (m.getIDCode() !== mol.getIDCode()) why ||= 'stereo changed';
    m.ensureHelperArrays(OCL.Molecule.cHelperCIP);
    const sides = []; for (let a = 0; a < mol.getAtoms(); a++) if ([1, 2].includes(m.getAtomCIPParity(a)) && !m.isAtomParityPseudo(a)) sides.push(Math.sign((P[a][0] - r.c[0]) * r.n[0] + (P[a][1] - r.c[1]) * r.n[1] + (P[a][2] - r.c[2]) * r.n[2]));
    if (!sides.includes(1) || !sides.includes(-1)) why ||= 'stereocenters not on both sides';
    ok &&= !why;
  }
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'meso plane', name.padEnd(46), got, r ? `(${r.how}, worst ${r.worst.toFixed(2)} Å)` : '', `${ms}ms`, ok ? '' : `(want ${want}${why ? `; ${why}` : ''})`);
}
{ const m = build3D(OCL, 'C[C@@H]1CCC[C@H](C)C1').mol; structuredClone(mirrorPlane(OCL, m)); } // the worker posts it

// ---- Fischer projection: [name, SMILES, rows top to bottom (left|right at each cross), D/L, meso look]
// A drawing turned 180° in the page is the same configuration, so that form also passes (meso-tartaric acid either way).
// R/S is read off the drawing itself (vertical groups back, horizontal ones forward, chem.js priorities) and must match
// chem.js's R/S from the 3D model, which test.mjs checks against RDKit above.
const FISCHER = [
  ['D-glyceraldehyde', 'C([C@H](C=O)O)O', 'CHO / H|OH / CH2OH', 'D'], ['L-alanine', 'C[C@@H](C(=O)O)N', 'COOH / NH2|H / CH3', 'L'],
  ['D-glucose', SMI['D-glucose'], 'CHO / H|OH / OH|H / H|OH / H|OH / CH2OH', 'D'], ['L-glucose', SMI['L-glucose'], 'CHO / OH|H / H|OH / OH|H / OH|H / CH2OH', 'L'],
  ['D-galactose', SMI['D-galactose'], 'CHO / H|OH / OH|H / OH|H / H|OH / CH2OH', 'D'], ['D-fructose', 'OC[C@@H](O)[C@@H](O)[C@H](O)C(=O)CO', 'CH2OH / C=O / OH|H / H|OH / H|OH / CH2OH', 'D'],
  ['L-serine', 'C([C@@H](C(=O)O)N)O', 'COOH / NH2|H / CH2OH', 'L'], ['L-threonine (L by C2, the α-carbon)', 'C[C@H]([C@@H](C(=O)O)N)O', 'COOH / NH2|H / H|OH / CH3', 'L'],
  ['L-alanine zwitterion', 'C[C@@H](C(=O)[O-])[NH3+]', 'COO⁻ / NH3⁺|H / CH3', 'L'],
  ['(2R,3R)-tartaric acid', SMI['(2R,3R)-tartaric acid'], 'COOH / H|OH / OH|H / COOH', ''], ['meso-tartaric acid', SMI['meso-tartaric acid'], 'COOH / H|OH / H|OH / COOH', '', true],
  ['meso-2,3-dibromobutane', 'C[C@@H](Br)[C@H](C)Br', 'CH3 / H|Br / H|Br / CH3', '', true], ['meso-2,4-pentanediol', 'C[C@@H](O)C[C@H](C)O', 'CH3 / H|OH / CH2 / H|OH / CH3', '', true],
  ['(R)-2-butanol', SMI['(R)-2-butanol'], 'CH3 / OH|H / CH2CH3', ''], ['(S)-2-butanol', SMI['(S)-2-butanol'], 'CH3 / H|OH / CH2CH3', ''],
  ['3-chlorobutan-2-ol: OH end on top', 'C[C@@H](O)[C@@H](C)Cl', 'CH3 / OH|H / H|Cl / CH3', ''], ['L-leucine', 'CC(C)C[C@@H](C(=O)O)N', 'COOH / NH2|H / CH2CH(CH3)2', 'L'],
  ['cis-1,2-dimethylcyclohexane', SMI['cis-1,2-dimethylcyclohexane'], 'ring', ''], ['menthol', 'C[C@@H]1CC[C@H]([C@@H](C1)O)C(C)C', 'ring', ''], ['ethanol', 'CCO', 'none', ''],
];
const fRows = f => f.why || f.rows.map(r => (r.cross ? `${r.left.label}|${r.right.label}` : r.label)).join(' / ');
const turned180 = s => s.split(' / ').reverse().map(x => (x.includes('|') ? x.split('|').reverse().join('|') : x)).join(' / ');
const IDEAL = { up: [0, 1, -1], down: [0, -1, -1], left: [-1, 0, 1], right: [1, 0, 1] };
for (const [name, smi, want, dl, meso = false] of FISCHER) {
  const { mol, given } = build3D(OCL, smi), rows = stereo(OCL, mol, given), f = fischer(OCL, mol, rows), got = fRows(f);
  const bad = [];
  if (got !== want && turned180(got) !== want) bad.push(`rows ${got}`);
  if ((f.dl?.dl || '') !== dl) bad.push(`D/L ${f.dl?.dl || 'none'}`);
  if (!f.why && (f.mirror !== null) !== meso) bad.push(`mirror ${f.mirror}`);
  for (const r of f.rows || []) {
    if (!r.cross || !r.rs) continue;
    const at = { [r.up]: IDEAL.up, [r.down]: IDEAL.down, [r.left.atom]: IDEAL.left, [r.right.atom]: IDEAL.right };
    const [v1, v2, v3, v4] = cipOrder(mol, r.atom).map(a => at[a]), s = (p, q) => p.map((x, k) => x - q[k]);
    const a = s(v1, v4), b = s(v2, v4), c = s(v3, v4), drawn = a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0]) < 0 ? 'R' : 'S';
    if (drawn !== r.rs) bad.push(`drawing reads ${drawn} at atom ${r.atom}, model ${r.rs}`);
  }
  const ok = !bad.length;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'fischer', name.padEnd(36), got, f.dl ? `[${f.dl.dl}]` : '', f.mirror !== null && f.mirror !== undefined ? '[mirror]' : '', ok ? '' : `(${bad.join('; ')}; want ${want} ${dl})`);
}
{ // the drawing comes from the configuration, so another conformation (D-glucose, seed 5) draws the same
  const b = build3D(OCL, SMI['D-glucose']), m = new OCL.ConformerGenerator(5).getOneConformerAsMolecule(OCL.Molecule.fromSmiles(SMI['D-glucose']));
  const one = fRows(fischer(OCL, b.mol, stereo(OCL, b.mol, b.given))), two = fRows(fischer(OCL, m, stereo(OCL, m, b.given)));
  if (one !== two) { fail++; console.log('FAIL fischer: D-glucose in another conformation draws', two, 'not', one); }
  structuredClone(fischer(OCL, b.mol, stereo(OCL, b.mol, b.given))); // the worker posts it
}

// ---- The same checks on the SMILES the site gets for these names (OPSIN, or PubChem where OPSIN fails), as returned on
// 17 Oct 2026: written differently from the ones above, so this also checks that nothing depends on how a SMILES is written.
const LIVE = {
  '(R)-2-butanol': 'C[C@H](CC)O', '(S)-2-butanol': 'C[C@@H](CC)O', '(2R,3R)-tartaric acid': '[C@@H]([C@H](C(=O)O)O)(C(=O)O)O',
  '(2S,3S)-tartaric acid': '[C@H]([C@@H](C(=O)O)O)(C(=O)O)O', 'meso-tartaric acid': '[C@@H]([C@@H](C(=O)O)O)(C(=O)O)O',
  'cis-2-butene': 'C\\C=C/C', 'trans-2-butene': 'C\\C=C\\C', 'D-glucose': 'O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO',
  'D-galactose': 'O=C[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO', 'L-glucose': 'O=C[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)CO',
  'cis-1,2-dimethylcyclohexane': 'C[C@H]1[C@H](CCCC1)C', 'trans-1,2-dimethylcyclohexane': 'C[C@H]1[C@@H](CCCC1)C',
  butane: 'CCCC', isobutane: 'CC(C)C', '1-propanol': 'C(CC)O', '2-propanol': 'CC(C)O', ethanol: 'C(C)O', 'dimethyl ether': 'COC', methanol: 'CO',
  'D-glyceraldehyde': 'O=C[C@H](O)CO', 'L-alanine': 'N[C@@H](C)C(=O)O', 'meso-2,3-dibromobutane': 'C[C@H]([C@H](C)Br)Br',
  'meso-2,4-pentanediol': 'C[C@H](C[C@H](C)O)O', 'cis-1,2-dimethylcyclopropane': 'C[C@H]1[C@H](C1)C', 'cis-1,3-dimethylcyclohexane': 'C[C@@H]1C[C@@H](CCC1)C',
};
for (const [x, y, want] of RELATION.filter(([x, y]) => LIVE[x] && LIVE[y])) {
  const got = relation(OCL, build3D(OCL, LIVE[x]), build3D(OCL, LIVE[y])).kind;
  if (got !== want) fail++;
  console.log(got === want ? 'ok  ' : 'FAIL', 'live relation', `${x} / ${y}`.padEnd(56), got, got === want ? '' : `(want ${want})`);
}
for (const [name, want] of [['meso-tartaric acid', true], ['meso-2,3-dibromobutane', true], ['meso-2,4-pentanediol', true], ['cis-1,2-dimethylcyclopropane', true],
  ['cis-1,3-dimethylcyclohexane', true], ['(2R,3R)-tartaric acid', false], ['(R)-2-butanol', false]]) {
  const got = !!mirrorPlane(OCL, build3D(OCL, LIVE[name]).mol);
  if (got !== want) fail++;
  console.log(got === want ? 'ok  ' : 'FAIL', 'live meso plane', name.padEnd(30), got ? 'plane' : 'no plane');
}
for (const [name, smi, want, dl] of FISCHER.filter(([n]) => LIVE[n])) {
  const { mol, given } = build3D(OCL, LIVE[name]), f = fischer(OCL, mol, stereo(OCL, mol, given)), got = fRows(f);
  const ok = (got === want || turned180(got) === want) && (f.dl?.dl || '') === dl;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'live fischer', name.padEnd(30), got, f.dl ? `[${f.dl.dl}]` : '', ok ? '' : `(want ${want} ${dl})`);
}

// ---- Shape tab: Newman bonds, one per symmetry class: [name, SMILES, 'front-back' bonds (* = anti/gauche references), ranked]
const NEWMAN = [
  ['ethane', 'CC', '0-1', true], ['propane', 'CCC', '0-1', true], ['butane', 'CCCC', '0-1 1-2*', true],
  ['2-methylbutane', 'CC(C)CC', '0-1 1-3 3-4', true], ['2,3-dimethylbutane', 'CC(C)C(C)C', '0-1 1-3', true],
  ['2,2-dimethylbutane', 'CC(C)(C)CC', '0-1 1-4 4-5', true], ['methylcyclohexane', 'CC1CCCCC1', '0-1', true],
  ['1,2-dichloroethane', 'ClCCCl', '1-2*', false], ['ethanol', 'CCO', '0-1', false],
  ['propene', 'CC=C', '', false], ['toluene', 'Cc1ccccc1', '', false], ['cyclohexane', 'C1CCCCC1', '', false],
  ['propylbenzene', 'CCCc1ccccc1', '0-1 1-2*', true], ['1-pentene', 'C=CCCC', '2-3* 3-4', true], ['propanal', 'CCC=O', '0-1', false],
  ['hexenylcyclohexanol', 'C(CC/C=C/C)[C@H]1[C@@H](CCCC1)O', '0-1* 1-2* 0-6', true],
];
// Independent check of the order: Newman angles straight from the 3D model (own projection, no rounding), each form's
// energy from the textbook increments, ranked the same way. The result must match `levels` exactly.
// Same chemistry rule as chem.js: any two non-H groups count when eclipsed, only two sp³ carbons when gauche.
const INC = { HH: 4.0, HX: 6.0, XX: 11, gauche: 3.8 }; // kJ/mol: eclipsed pairs, and gauche CH3/CH3
const vsub = (p, q) => p.map((v, i) => v - q[i]), vdot = (p, q) => p.reduce((s, v, i) => s + v * q[i], 0);
function textbookLevels(mol, at, { atoms: [F, B], front, back, x, phi0 }) {
  const P = a => [mol.getAtomX(a), mol.getAtomY(a), mol.getAtomZ(a)];
  const ax = vsub(P(B), P(F)), u = ax.map(v => v / Math.hypot(...ax));
  const flat = a => { const v = vsub(P(a), P(front.includes(a) ? F : B)); return vsub(v, u.map(k => k * vdot(v, u))); };
  const r = flat(x[0]), s = [u[1] * r[2] - u[2] * r[1], u[2] * r[0] - u[0] * r[2], u[0] * r[1] - u[1] * r[0]];
  const th = a => { const v = flat(a); return Math.atan2(vdot(v, s), vdot(v, r)) * 180 / Math.PI; };
  const heavy = a => mol.getAtomicNo(a) !== 1, sp3 = a => at[a].el === 'C' && at[a].hyb === 'sp3';
  const close = []; // pairs in steric strain at each form: two sp³ carbons eclipsed or gauche, any two non-H groups eclipsed
  const E = [0, 60, 120, 180, 240, 300].map(phi => {
    let e = 0, n = 0;
    for (const f of front) for (const b of back) {
      const d = Math.abs(((th(b) + phi - phi0 - th(f)) % 360 + 540) % 360 - 180); // 0 = eclipsed, 180 = anti
      if ((d < 90 && sp3(f) && sp3(b)) || (d < 30 && heavy(f) && heavy(b))) n++;
      if (d < 30) e += heavy(f) && heavy(b) ? INC.XX : heavy(f) || heavy(b) ? INC.HX : INC.HH;
      else if (d < 90 && sp3(f) && sp3(b)) e += INC.gauche;
    }
    close.push(n);
    return Math.round(e * 10);
  });
  const keys = [...new Set(E)].sort((p, q) => p - q);
  return Object.assign(E.map(k => keys.indexOf(k)), { close });
}
for (const [name, smi, want, ranked] of NEWMAN) {
  const { mol } = build3D(OCL, smi), at = atomInfo(mol), bonds = newmanBonds(OCL, mol, at);
  const got = bonds.map(b => b.atoms.join('-') + (b.ref ? '*' : '')).join(' ');
  // ranked: order and steric pairs at each form both match the independent check; unranked: no order, no steric pairs
  const bad = bonds.filter(b => ranked
    ? b.levels?.join() !== textbookLevels(mol, at, b).join() || [0, 60, 120, 180, 240, 300].map(p => stericPairs(b, p).length).join() !== textbookLevels(mol, at, b).close.join()
    : b.levels || [0, 30, 60].some(p => stericPairs(b, p).length));
  const ok = got === want && !bad.length;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'newman', name.padEnd(20), got || '(none)', bonds.map(b => b.levels?.join('') ?? '-').join(' '),
    ok ? '' : `(want ${want}, ${ranked ? 'textbook ' + bad.map(b => textbookLevels(mol, at, b).join('')).join(' ') : 'no ranking'})`);
}
// Butane C2–C3 by hand: 0° (CH3/CH3 eclipsed) > 120° = 240° (H/CH3 eclipsed) > 60° = 300° (gauche) > 180° (anti)
const bu = newmanBonds(OCL, build3D(OCL, 'CCCC').mol, atomInfo(build3D(OCL, 'CCCC').mol))[1];
if (bu.levels.join('') !== '312021' || Math.abs(bu.phi0 - 180) > 2) { fail++; console.log('FAIL butane C2–C3 order', bu.levels, bu.phi0); }
// Butane C2–C3 steric pairs by slider angle: the two CH3 carbons (atoms 0 and 3) from gauche to eclipsed, none from 70° to 290°
for (const [phi, want] of [[0, 1], [30, 1], [60, 1], [65, 1], [70, 0], [120, 0], [180, 0], [240, 0], [295, 1], [300, 1]]) {
  const got = stericPairs(bu, phi);
  if (got.length !== want || (want && got[0].slice(0, 2).sort().join() !== '0,3')) { fail++; console.log('FAIL butane steric pairs at', phi, got, `(want ${want})`); }
}
const ethB = newmanBonds(OCL, build3D(OCL, 'CC').mol, atomInfo(build3D(OCL, 'CC').mol))[0];
if ([0, 60].some(p => stericPairs(ethB, p).length)) { fail++; console.log('FAIL ethane has no steric pairs'); }

// Propylbenzene CH2–CH2 (phenyl and CH3): eclipsed groups still highest, but anti and gauche at one level, and the
// pink rod only when eclipsed. Pairs: the CH3 carbon (atom 0) and the ring carbon (atom 3).
const pb = newmanBonds(OCL, build3D(OCL, 'CCCc1ccccc1').mol, atomInfo(build3D(OCL, 'CCCc1ccccc1').mol))[1];
if (pb.levels.join('') !== '201010' || !pb.flat) { fail++; console.log('FAIL propylbenzene order', pb.levels, pb.flat); }
for (const [phi, want] of [[0, 1], [5, 1], [10, 0], [60, 0], [180, 0]])
  if (stericPairs(pb, phi).length !== want) { fail++; console.log('FAIL propylbenzene steric pairs at', phi, `(want ${want})`); }

// Why a flat group counts when eclipsed but not when gauche: MMFF94s+ minima from anti and gauche starts.
// gauche − anti for butane (CH3/CH3) is about 3.3 kJ/mol; with a phenyl or vinyl group it is under 1 kJ/mol.
function gaucheMinusAnti(smi, [a, b, c, d]) {
  const { mol } = build3D(OCL, smi);
  const P = (m, i) => [m.getAtomX(i), m.getAtomY(i), m.getAtomZ(i)];
  const crs = (p, q) => [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]];
  const dih = m => { const [pa, pb, pc, pd] = [a, b, c, d].map(i => P(m, i)), b1 = vsub(pb, pa), b2 = vsub(pc, pb), b3 = vsub(pd, pc);
    const n1 = crs(b1, b2), n2 = crs(b2, b3); return Math.atan2(vdot(crs(n1, b2.map(v => v / Math.hypot(...b2))), n2), vdot(n1, n2)) * 180 / Math.PI; };
  const energyAt = target => {
    const m = mol.getCompactCopy(); m.ensureHelperArrays(OCL.Molecule.cHelperNeighbours);
    const side = new Set([c]); for (const i of side) for (let k = 0; k < m.getAllConnAtoms(i); k++) { const j = m.getConnAtom(i, k); if (j !== b) side.add(j); }
    const o = P(m, b), ax = vsub(P(m, c), o), u = ax.map(v => v / Math.hypot(...ax)), t = (target - dih(m)) * Math.PI / 180;
    for (const i of side) { const v = vsub(P(m, i), o), kv = crs(u, v), kd = vdot(u, v);
      const r = v.map((x, j) => x * Math.cos(t) + kv[j] * Math.sin(t) + u[j] * kd * (1 - Math.cos(t)) + o[j]); m.setAtomX(i, r[0]); m.setAtomY(i, r[1]); m.setAtomZ(i, r[2]); }
    new OCL.ForceFieldMMFF94(m, 'MMFF94s+').minimise();
    return new OCL.ForceFieldMMFF94(m, 'MMFF94s+').getTotalEnergy() * 4.184;
  };
  return Math.min(energyAt(60), energyAt(-60)) - energyAt(180);
}
for (const [name, smi, quad, lo, hi] of [['butane', 'CCCC', [0, 1, 2, 3], 2, 5], ['propylbenzene', 'CCCc1ccccc1', [0, 1, 2, 3], -0.5, 1.5],
  ['1-pentene', 'C=CCCC', [1, 2, 3, 4], -0.5, 1.5]]) {
  const d = gaucheMinusAnti(smi, quad), ok = d >= lo && d <= hi;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'MMFF gauche − anti', name.padEnd(14), `${d.toFixed(2)} kJ/mol`, ok ? '' : `(want ${lo} to ${hi})`);
}

// What the slider position shows: [phi, references?, expected]
for (const [phi, ref, want] of [[180, true, 'staggered anti'], [60, true, 'staggered gauche'], [300, true, 'staggered gauche'],
  [0, true, 'eclipsed'], [355, true, 'eclipsed'], [120, false, 'eclipsed'], [63, true, 'staggered gauche'], [66, true, 'between'],
  [30, false, 'between'], [180, false, 'staggered']]) {
  const c = conformer(phi, ref), got = c.kind ? [c.kind, c.rel].filter(Boolean).join(' ') : 'between';
  if (got !== want) { fail++; console.log('FAIL conformer at', phi, ref, got, `(want ${want})`); }
}

// Ring angles: [name, SMILES, ring size, lowest and highest allowed mean angle, ideal]
const RINGS = [
  ['cyclopropane', 'C1CC1', 3, 59.5, 60.5, 109.5], ['cyclobutane', 'C1CCC1', 4, 86, 90, 109.5],
  ['cyclopentane', 'C1CCCC1', 5, 103, 106, 109.5], ['cyclohexane', 'C1CCCCC1', 6, 110, 112.5, 109.5],
  ['benzene', 'c1ccccc1', 6, 119.5, 120.5, 120], ['oxirane', 'C1CO1', 3, 59, 61, 109.5],
];
for (const [name, smi, size, lo, hi, ideal] of RINGS) {
  const { mol } = build3D(OCL, smi), r = ringAngles(mol, atomInfo(mol));
  const ok = r.length === 1 && r[0].size === size && r[0].mean >= lo && r[0].mean <= hi && Math.abs(r[0].ideal - ideal) < 0.01;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'ring', name.padEnd(14), r.map(x => `${x.size}: mean ${x.mean.toFixed(1)}°, ideal ${x.ideal}`).join(' '), ok ? '' : `(want ${size}: ${lo}–${hi}°, ideal ${ideal})`);
}
const decalin = ringAngles(build3D(OCL, 'C1CCC2CCCCC2C1').mol, atomInfo(build3D(OCL, 'C1CCC2CCCCC2C1').mol));
if (decalin.length !== 2 || decalin.some(r => r.size !== 6)) { fail++; console.log('FAIL decalin rings', decalin.map(r => r.size)); }

// ---- Chairs: [name, SMILES, flip block of each chair found ('-' = can flip), every non-H group equatorial?]
const GLUCOSE = 'C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O', MENTHOL = 'C[C@@H]1CC[C@H]([C@@H](C1)O)C(C)C';
const CHAIRS = [
  ['cyclohexane', 'C1CCCCC1', 'h'], ['piperidine', 'C1CCNCC1', 'h'], ['tetrahydropyran', 'C1CCOCC1', 'h'], ['1,4-dioxane', 'C1COCCO1', 'h'],
  ['methylcyclohexane', 'CC1CCCCC1', '-'], ['N-methylpiperidine', 'CN1CCCCC1', '-'],
  ['cyclohexanone', 'O=C1CCCCC1', ''], ['cyclohexene', 'C1CCC=CC1', ''], ['benzene', 'c1ccccc1', ''], ['cyclopentane', 'C1CCCC1', ''],
  ['cycloheptane', 'C1CCCCCC1', ''], ['norbornane (boat)', 'C1CC2CCC1C2', ''], ['bicyclo[2.2.2]octane (boats)', 'C1CC2CCC1CC2', ''],
  ['trans-decalin', 'C1CC[C@H]2CCCC[C@@H]2C1', 'trans trans'], ['cis-decalin', 'C1CC[C@H]2CCCC[C@H]2C1', 'cis'],
  ['adamantane', 'C1C2CC3CC1CC(C2)C3', 'bridged bridged bridged bridged'], ['spiro[5.5]undecane', 'C1CCC2(CC1)CCCCC2', '- -'],
  ['beta-D-glucose (all equatorial)', GLUCOSE, '-', true], ['(1R,2S,5R)-menthol (all equatorial)', MENTHOL, '-', true],
];
// Geometry every chair must have: each ring atom with two outside atoms has one axial and one equatorial, on opposite
// faces; axial bonds alternate up and down around the ring.
const chairOk = c => c.atoms.every((a, i) => {
  const s = c.subs.filter(x => x.atom === a), nx = c.subs.find(x => x.atom === c.atoms[(i + 1) % 6] && x.ax), ax = s.find(x => x.ax);
  return (s.length !== 2 || (s.filter(x => x.ax).length === 1 && s[0].up !== s[1].up)) && (!ax || !nx || ax.up !== nx.up);
});
for (const [name, smi, want, allEq] of CHAIRS) {
  const { mol } = build3D(OCL, smi), ch = chairs(mol, atomInfo(mol));
  const got = ch.map(c => c.block || '-').join(' ');
  const ok = got === want && ch.every(chairOk) && (!allEq || ch.every(c => c.subs.every(s => s.h || !s.ax)));
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'chair', name.padEnd(34), got || '(none)', ok ? '' : `(want ${want || 'none'}${allEq ? ', all equatorial' : ''}${ch.every(chairOk) ? '' : ', bad axial/equatorial geometry'})`);
}

// Axial/equatorial against the textbook table, using the ring cis/trans that stereo() reads (checked against RDKit above):
// cis-1,2, trans-1,3 and cis-1,4 put one group axial and one equatorial; the others put both on the same kind.
// Also cis = the two groups on one face.
const DIMETHYL = [
  ['cis-1,2', 'C[C@H]1CCCC[C@H]1C', 'cis', 'ae'], ['trans-1,2', 'C[C@H]1CCCC[C@@H]1C', 'trans', 'same'],
  ['cis-1,3', 'C[C@@H]1CCC[C@H](C)C1', 'cis', 'same'], ['trans-1,3', 'C[C@@H]1CCC[C@@H](C)C1', 'trans', 'ae'],
  ['cis-1,4', 'C[C@H]1CC[C@@H](C)CC1', 'cis', 'ae'], ['trans-1,4', 'C[C@H]1CC[C@H](C)CC1', 'trans', 'same'],
];
const axeq = subs => { const g = subs.filter(s => !s.h); return { kind: g[0].ax === g[1].ax ? 'same' : 'ae', face: g[0].up === g[1].up ? 'cis' : 'trans', g }; };
for (const [name, smi, ct, want] of DIMETHYL) {
  const { mol, given } = build3D(OCL, smi), [c] = chairs(mol, atomInfo(mol)), ringCT = stereo(OCL, mol, given).find(r => r.kind === 'ring')?.ct;
  const got = axeq(c.subs), ok = ringCT === ct && got.kind === want && got.face === ct && chairOk(c) && !c.block && !c.judge;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'ax/eq', `${name}-dimethyl`.padEnd(20), got.g.map(s => (s.ax ? 'axial' : 'equatorial')).join(', '), ok ? '' : `(want ${ct}: ${want}; stereo says ${ringCT}, faces ${got.face})`);
}

// 1,3-diaxial pairs, checked against distance: for an axial group X, the partner on each ring atom two away must be that
// atom's outside atom nearest to X. [SMILES, pairs expected when X is axial in the model]
const dist = (m, a, b) => Math.hypot(m.getAtomX(a) - m.getAtomX(b), m.getAtomY(a) - m.getAtomY(b), m.getAtomZ(a) - m.getAtomZ(b));
for (const [smi, want] of [['CC1CCCCC1', 2], ['OC1CCCCC1', 2], ['C[C@@H]1CCC[C@H](C)C1', 3], ['CC1(C)CCCCC1', 2], ['C1CCCCC1', 0], [MENTHOL, 0]]) {
  const { mol } = build3D(OCL, smi), [c] = chairs(mol, atomInfo(mol)), pairs = diaxial(c);
  const axial = c.subs.filter(s => s.ax && !s.h), n = axial.length ? want : 0;
  const near = axial.every(s => { const i = c.atoms.indexOf(s.atom);
    return [2, 4].every(d => { const r = c.atoms[(i + d) % 6], o = c.subs.filter(x => x.atom === r).sort((p, q) => dist(mol, p.sub, s.sub) - dist(mol, q.sub, s.sub))[0];
      return pairs.some(([x, y]) => (x === s.sub && y === o.sub) || (y === s.sub && x === o.sub)); }); });
  const ok = pairs.length === n && near;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', '1,3-diaxial', smi.padEnd(30), axial.length ? `${axial.length} axial group(s), ${pairs.length} pairs` : 'no axial group', ok ? '' : `(want ${n} pairs, nearest-atom check ${near})`);
}

// ---- Ring flip: the other chair from OCL conformers + MMFF94s+. [name, SMILES, flip expected?]
// tert-butyl: OCL's generator never builds the axial chair (its ring-fragment step keeps only shapes close in energy),
// so the page says "flip not available"; this row checks that it stays honest about it.
const FLIPS = [['methylcyclohexane', 'CC1CCCCC1', true], ['cyclohexanol', 'OC1CCCCC1', true], ['chlorocyclohexane', 'ClC1CCCCC1', true],
  ['cyclohexylamine', 'NC1CCCCC1', true], ['1,1-dimethylcyclohexane', 'CC1(C)CCCCC1', true], ['N-methylpiperidine', 'CN1CCCCC1', true],
  ['fluorocyclohexane', 'FC1CCCCC1', true], ['bromocyclohexane', 'BrC1CCCCC1', true], ['iodocyclohexane', 'IC1CCCCC1', true], ['cyclohexanecarbonitrile', 'N#CC1CCCCC1', true],
  ...DIMETHYL.map(([n, s]) => [`${n}-dimethyl`, s, true]), ['tert-butylcyclohexane', 'CC(C)(C)C1CCCCC1', false]];
const flipped = {}; // name → [flip result, model chair], reused by the stability checks below
for (const [name, smi, want] of FLIPS) {
  const { mol } = build3D(OCL, smi), [c] = chairs(mol, atomInfo(mol)), t = performance.now(), f = flipChair(OCL, mol, c.atoms);
  flipped[name] = [f, c];
  let why = '';
  if (!!f !== want) why = want ? 'no flipped chair found' : 'found a flip';
  else if (f) {
    // the shown coordinates: an MMFF94s+ minimum (minimise() stays put), moved rigidly onto the model's ring
    const m = mol.getCompactCopy(); m.ensureHelperArrays(OCL.Molecule.cHelperNeighbours);
    f.xyz.forEach(([x, y, z], i) => { m.setAtomX(i, x); m.setAtomY(i, y); m.setAtomZ(i, z); });
    const ff = new OCL.ForceFieldMMFF94(m, 'MMFF94s+'), e0 = ff.getTotalEnergy();
    const P = (mm, a) => [mm.getAtomX(a), mm.getAtomY(a), mm.getAtomZ(a)], cen = mm => [0, 1, 2].map(k => c.atoms.reduce((s, a) => s + P(mm, a)[k] / 6, 0));
    const [g0, g1] = [c.subs, chairs(m, atomInfo(m))[0]?.subs ?? []].map(s => s.filter(x => !x.h));
    if (ff.minimise() !== 0 || Math.abs(ff.getTotalEnergy() - e0) > 0.01) why = 'shown coordinates are not an MMFF minimum';
    else if (Math.hypot(...vsub(cen(m), cen(mol))) > 0.01) why = 'ring not moved onto the model ring';
    else if (g1.length !== g0.length || !g1.every((s, i) => s.ax !== g0[i].ax && s.up === g0[i].up)) why = 'groups not swapped axial/equatorial on the same faces';
    else if (!f.subs.filter(x => !x.h).every((s, i) => s.ax === g1[i].ax)) why = 'returned axial/equatorial does not match the coordinates';
    else if (!c.subs.filter(x => x.h).every(s => chairs(m, atomInfo(m))[0].subs.find(x => x.sub === s.sub).up === s.up)) why = 'a ring H moved to the other face (H identity lost)';
  }
  if (why) fail++;
  console.log(why ? 'FAIL' : 'ok  ', 'flip', name.padEnd(26), f ? f.subs.filter(x => !x.h).map(s => (s.ax ? 'axial' : 'equatorial')).join(', ') : 'not available',
    `${(performance.now() - t).toFixed(0)}ms`, why ? `(${why})` : '');
}

// ---- More stable chair. Textbook: a lone group prefers equatorial; cis-1,2, trans-1,3 and cis-1,4 have one axial group
// in either chair (same energy); trans-1,2, cis-1,3 and trans-1,4 prefer both equatorial.
const STABLE = [['methylcyclohexane', 'eq'], ['cyclohexanol', 'eq'], ['cyclohexylamine', 'eq'], ['N-methylpiperidine', 'eq'],
  ['cis-1,2-dimethyl', 'same'], ['trans-1,2-dimethyl', 'eq'], ['cis-1,3-dimethyl', 'eq'], ['trans-1,3-dimethyl', 'same'],
  ['cis-1,4-dimethyl', 'same'], ['trans-1,4-dimethyl', 'eq'], ['1,1-dimethylcyclohexane', 'same']];
for (const [name, want] of STABLE) {
  const [f, c] = flipped[name], v = moreStable(f.dE), win = v === 'built' ? c.subs : f.subs;
  const got = v === 'same' ? 'same' : win.filter(s => !s.h).every(s => !s.ax) ? 'eq' : 'not all equatorial';
  const ok = got === want && !c.judge;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'stable', name.padEnd(26), v === 'same' ? 'both chairs about the same' : 'all-equatorial chair lower', `(ΔE ${Math.abs(f.dE).toFixed(2)} kJ/mol)`, ok ? '' : `(want ${want}${c.judge ? ', but judged ' + c.judge.why : ''})`);
}
// MMFF94s+ Δ(axial − equatorial) for one group against measured A-values (kJ/mol; Wikipedia's table in kcal/mol × 4.184,
// with the brief's values for CH₃, OH, Cl). Groups the page judges must come out equatorial-favoured and in the A-value
// order (CH₃ > NH₂ > OH); groups it leaves out must still disagree, so a fixed force field would show up here.
// tert-butyl (A ≈ 20) cannot be checked: no axial chair from OCL.
const A = { methylcyclohexane: 7.3, cyclohexylamine: 6.7, cyclohexanol: 3.6, chlorocyclohexane: 2.2, fluorocyclohexane: 0.63, bromocyclohexane: 1.6, iodocyclohexane: 1.8, cyclohexanecarbonitrile: 0.71 };
const axMinusEq = name => { const [f, c] = flipped[name]; return c.subs.find(s => !s.h).ax ? -f.dE : f.dE; };
for (const name of Object.keys(A)) {
  const d = axMinusEq(name), judged = !flipped[name][1].judge, eqSide = d >= 1;
  const ok = judged ? eqSide : ['halogen', 'nitrile'].includes(flipped[name][1].judge.why) && (!eqSide || name === 'iodocyclohexane'); // I: right side, left out with the other halogens
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'A-value', name.padEnd(24), `MMFF ${d.toFixed(2)} vs A ${A[name]} kJ/mol`, judged ? 'judged' : 'left out (MMFF disagrees)', ok ? '' : '(recheck noVerdict in chem.js)');
}
if (!(axMinusEq('methylcyclohexane') > axMinusEq('cyclohexylamine') && axMinusEq('cyclohexylamine') > axMinusEq('cyclohexanol'))) { fail++; console.log('FAIL MMFF order CH₃ > NH₂ > OH'); }

// ---- Polarity: textbook polar / nonpolar from the vector sum of bond dipoles (POLAR_BOND, POLAR_NET in chem.js).
// Nonpolar ones must sum to under 0.02 and polar ones to 0.3 or more, so the 0.1 cutoff sits well clear of both.
// hexane has no polar bond at all, so its answer does not depend on the conformer.
const POL = [
  ['CO2', 'O=C=O', false], ['BF3', 'FB(F)F', false], ['CH4', 'C', false], ['CCl4', 'ClC(Cl)(Cl)Cl', false],
  ['trans-1,2-dichloroethene', 'Cl/C=C/Cl', false], ['p-dichlorobenzene', 'Clc1ccc(Cl)cc1', false], ['benzene', 'c1ccccc1', false],
  ['hexane', 'CCCCCC', false], ['CF4', 'FC(F)(F)F', false], ['1,3,5-trichlorobenzene', 'Clc1cc(Cl)cc(Cl)c1', false],
  ['H2O', 'O', true], ['NH3', 'N', true], ['CH3Cl', 'CCl', true], ['CH2Cl2', 'ClCCl', true], ['CHCl3', 'ClC(Cl)Cl', true],
  ['cis-1,2-dichloroethene', 'Cl/C=C\\Cl', true], ['o-dichlorobenzene', 'Clc1ccccc1Cl', true], ['acetone', 'CC(C)=O', true],
  ['ethanol', 'CCO', true], ['pyridine', 'c1ccncc1', true], ['acetonitrile', 'CC#N', true], ['trimethylamine', 'CN(C)C', true],
  ['bromomethane', 'CBr', true],
  ['iodomethane', 'CI', null], ['dimethyl sulfide', 'CSC', null], ['acetate', 'CC(=O)[O-]', null],
];
for (const [name, smi, want] of POL) {
  const { mol, minimized } = build3D(OCL, smi), p = polarity(mol), len = Math.hypot(...p.net);
  // BF3: MMFF94s+ has no boron types, so the raw conformer (angles 118–121°) stays and sums to about 0.06
  const ok = p.polar === want && (want === null || (want ? len >= 0.3 : len < (minimized ? 0.02 : 0.08))) && (name !== 'hexane' || !p.bonds.length);
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'polar', name.padEnd(26), p.polar === null ? `no verdict (${p.why})` : p.polar ? 'polar' : 'nonpolar', `sum ${len.toFixed(3)}`, ok ? '' : `(want ${want})`);
}
// Conformation decides for 1,2-dichloroethane: anti (Cl–C–C–Cl 180°) sums to 0, gauche does not.
{ const { mol } = build3D(OCL, 'ClCCCl'), p = polarity(mol);
  const [a, b] = [0, 3], t = Math.abs(dihedralOf(mol, a, 1, 2, b));
  const ok = t > 170 && p.polar === false;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'polar 1,2-dichloroethane anti in the model:', `Cl–C–C–Cl ${t.toFixed(0)}°`, p.polar ? 'polar' : 'nonpolar');
  // turn one CH2Cl by 120° about the C–C bond (to gauche) and sum again
  const g = mol.getCompactCopy(), [o, u] = [[g.getAtomX(1), g.getAtomY(1), g.getAtomZ(1)], [g.getAtomX(2), g.getAtomY(2), g.getAtomZ(2)]];
  const k = u.map((x, i) => x - o[i]), n = Math.hypot(...k), ax = k.map(x => x / n), c = Math.cos(2 * Math.PI / 3), s = Math.sin(2 * Math.PI / 3);
  g.ensureHelperArrays(OCL.Molecule.cHelperNeighbours);
  const side = [3, ...Array.from({ length: g.getAllConnAtoms(2) }, (_, j) => g.getConnAtom(2, j)).filter(x => x !== 1)];
  for (const i of new Set(side)) {
    const v = [g.getAtomX(i) - u[0], g.getAtomY(i) - u[1], g.getAtomZ(i) - u[2]], d = v[0] * ax[0] + v[1] * ax[1] + v[2] * ax[2];
    const cr = [ax[1] * v[2] - ax[2] * v[1], ax[2] * v[0] - ax[0] * v[2], ax[0] * v[1] - ax[1] * v[0]];
    const r = v.map((x, j) => x * c + cr[j] * s + ax[j] * d * (1 - c));
    g.setAtomX(i, r[0] + u[0]); g.setAtomY(i, r[1] + u[1]); g.setAtomZ(i, r[2] + u[2]);
  }
  const pg = polarity(g), t2 = Math.abs(dihedralOf(g, a, 1, 2, b)), ok2 = Math.abs(t2 - 60) < 10 && pg.polar === true;
  if (!ok2) fail++;
  console.log(ok2 ? 'ok  ' : 'FAIL', 'polar 1,2-dichloroethane gauche:', `Cl–C–C–Cl ${t2.toFixed(0)}°`, pg.polar ? 'polar' : 'nonpolar');
}

// ---- Aromatic rings: π electrons per ring, per atom (heavy atoms in SMILES order), and per fused system.
// [name, SMILES, rings as "π: atom contributions" sorted, system totals]
const AROM = [
  ['benzene', 'c1ccccc1', ['6:111111'], [6]],
  ['pyridine', 'c1ccncc1', ['6:111111'], [6]], // N: C=N bond 1, its lone pair in the plane
  ['pyrrole', 'c1cc[nH]c1', ['6:11121'], [6]],
  ['furan', 'c1ccoc1', ['6:11121'], [6]],
  ['thiophene', 'c1ccsc1', ['6:11121'], [6]],
  ['imidazole', 'c1cnc[nH]1', ['6:11112'], [6]],
  ['naphthalene', 'c1ccc2ccccc2c1', ['6:111111', '6:111111'], [10]],
  ['indole', 'c1ccc2[nH]ccc2c1', ['6:111111', '6:12111'], [10]], // 5-ring atoms 3 (C), 4 (N–H), 5, 6, 7
  ['cyclopentadienyl anion', '[cH-]1cccc1', ['6:21111'], [6]],
  ['tropylium', '[cH+]1cccccc1', ['6:0111111'], [6]],
];
for (const [name, smi, rings, systems] of AROM) {
  const { mol } = build3D(OCL, smi), at = atomInfo(mol), ar = aromaticRings(mol);
  const got = ar.rings.map(r => `${r.pi}:${[...r.parts].sort((x, y) => x.atom - y.atom).map(p => p.e).join('')}`).sort();
  const sys = ar.systems.map(s => s.pi).sort((x, y) => x - y);
  // 1a agreement: an atom giving a lone pair (2) is sp² by resonance; one with a double bond (1) is sp² without it
  const clash = ar.rings.flatMap(r => r.parts).filter(p => at[p.atom].hyb && (at[p.atom].hyb !== 'sp2' || !!at[p.atom].resonance !== (p.why === 'lp')));
  const planarN = name === 'pyridine' ? ar.rings[0].parts.find(p => at[p.atom].el === 'N').planar : true;
  const ok = got.join() === [...rings].sort().join() && sys.join() === systems.join() && !clash.length && planarN
    && ar.rings.every(r => huckel(ar.systems[r.system].pi) !== null);
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'arom', name.padEnd(24), got.join(' '), '| system', sys.join(' '), ok ? '' : `(want ${rings.join(' ')} | ${systems.join(' ')}${clash.length ? '; disagrees with hybridization at ' + clash.map(p => p.atom) : ''})`);
}
if (aromaticRings(build3D(OCL, 'C1=CC=CCC1').mol).rings.length) { fail++; console.log('FAIL arom: 1,3-cyclohexadiene is not aromatic'); }

// ---- Equivalent H (¹H NMR signals)
const NMR = [
  ['ethanol', 'CCO', 3], ['acetone', 'CC(C)=O', 1], ['ethyl acetate', 'CCOC(C)=O', 3], ['toluene', 'Cc1ccccc1', 4],
  ['p-xylene', 'Cc1ccc(C)cc1', 2], ['2-butanone', 'CCC(C)=O', 3], ['propanoic acid', 'CCC(=O)O', 3],
  ['2-chlorobutane', 'CCC(C)Cl', 5], ['(R)-2-chlorobutane', 'CC[C@@H](C)Cl', 5], ['cyclohexane', 'C1CCCCC1', 1],
  ['benzene', 'c1ccccc1', 1], ['isopropyl alcohol', 'CC(C)O', 3], ['1-chloropropane', 'CCCCl', 3],
  ['propene', 'C=CC', 4], ['methylcyclohexane', 'CC1CCCCC1', 8], ['meso-2,3-dibromobutane', 'C[C@@H](Br)[C@H](C)Br', 2],
];
for (const [name, smi, want] of NMR) {
  const { mol } = build3D(OCL, smi), { h, cls } = equivalentH(mol), n = new Set(cls).size;
  let ok = n === want && h.length === cls.length;
  if (name === '2-chlorobutane') { // the CH₂ H pair (on heavy atom 1) must be split: diastereotopic
    const pair = h.map((a, i) => [a, cls[i]]).filter(([a]) => mol.getConnAtom(a, 0) === 1).map(([, c]) => c);
    ok &&= pair.length === 2 && pair[0] !== pair[1];
  }
  if (name === 'ethanol') { // the CH₂ pair is enantiotopic: one signal
    const pair = h.map((a, i) => [a, cls[i]]).filter(([a]) => mol.getConnAtom(a, 0) === 1).map(([, c]) => c);
    ok &&= pair.length === 2 && pair[0] === pair[1];
  }
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'nmr', name.padEnd(24), `${n} signals`, ok ? '' : `(want ${want})`);
}

// ---- Reactions tab. Products are compared with the textbook product as RDKit canonical isomeric SMILES; R/S and E/Z
// come from RDKit's CIP labels on the product SMILES that chem.js returns.
const rdCIP = s => { const m = RDKit.get_mol(s), t = JSON.parse(m.get_stereo_tags()); m.delete(); return t; };
// A product as RDKit reads it from the molfile chem.js posts: stereo from the 3D coordinates, H removed, canonical SMILES.
const rdMol = mf => { const m = RDKit.get_mol(mf, JSON.stringify({ removeHs: true })), c = m.get_smiles(); m.delete(); return c; };
const rxMol = smi => { const b = build3D(OCL, smi); return { ...b, atoms: atomInfo(b.mol) }; };
const rxHal = smi => { const b = rxMol(smi); return { ...b, hal: halides(OCL, b.mol, b.atoms) }; };
// SN2 backside cone: more atoms in the nucleophile's path as the carbon gets more crowded
{ const n = smi => rxHal(smi).hal.rows[0].cone.length;
  const [me, et, ipr, tbu, neo] = ['CBr', 'CCBr', 'CC(C)Br', 'CC(C)(C)Br', 'CC(C)(C)CBr'].map(n);
  const ok = me < et && et < ipr && ipr < tbu && neo > et;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'SN2 cone atoms', `bromomethane ${me} < bromoethane ${et} < 2-bromopropane ${ipr} < 2-bromo-2-methylpropane ${tbu}; 1-bromo-2,2-dimethylpropane ${neo} > ${et}`); }
// [name, SMILES (the C–X center written [C@H](Br) or [C@@H](Br)), rate, product CIP at the C–O carbon ('' = not a stereocenter)]
// Every product must equal the reactant SMILES with Br → O and that center's @/@@ inverted. 3° and neopentyl: no product shown.
const SN2 = [
  ['bromomethane', 'CBr', 'fast', ''], ['bromoethane', 'CCBr', 'fast', ''], ['1-bromopropane', 'CCCBr', 'fast', ''], ['2-bromopropane', 'CC(C)Br', 'slow', ''],
  ['2-bromo-2-methylpropane', 'CC(C)(C)Br', 'none', ''], ['1-bromo-2,2-dimethylpropane', 'CC(C)(C)CBr', 'none', ''],
  ['(R)-2-bromobutane → (S)-2-butanol', 'C[C@@H](Br)CC', 'slow', 'S'], ['(S)-2-bromobutane → (R)-2-butanol', 'C[C@H](Br)CC', 'slow', 'R'],
  ['(R)-1-bromo-1-chloroethane → (R)-1-chloroethanol', 'C[C@@H](Br)Cl', null, 'R'], ['bromocyclohexane', 'BrC1CCCCC1', 'slow', ''],
  ['cis-1-bromo-2-methylcyclohexane → trans-2-methylcyclohexanol', 'C[C@H]1CCCC[C@@H](Br)1', 'slow', 'S'],
];
const invertBr = s => s.replace(/\[C(@@?)H\](?=\(Br\)|Br)/, m => m.replace(/@@?/, x => (x === '@' ? '@@' : '@'))).replace('Br', 'O');
for (const [name, smi, rate, cip] of SN2) {
  const { mol, hal, given } = rxHal(smi), h = hal.rows.find(x => x.el === 'Br'), r = sn2(OCL, mol, h.c, h.x), bad = [];
  if (hal.rows.length !== (smi.includes('Cl') ? 2 : 1)) bad.push(`${hal.rows.length} C–X rows`);
  if (sn2Rate(h) !== rate) bad.push(`rate ${sn2Rate(h)}`);
  if (!r) bad.push('no product');
  else {
    if (rdMol(r.product.molfile) !== rdCan(invertBr(smi))) bad.push(`product ${rdMol(r.product.molfile)}, want ${invertBr(smi)}`);
    const t = rdCIP(invertBr(smi)).CIP_atoms.find(([i]) => i === h.c)?.[1][1] ?? ''; // heavy atoms keep the SMILES order
    if (t !== cip || (r.product.rs ?? '') !== cip) bad.push(`R/S ${r.product.rs} (RDKit ${t || 'none'})`);
  }
  if (name.startsWith('cis-1-bromo')) { // the ring pair turns from cis to trans
    const was = stereo(OCL, mol, given).find(x => x.kind === 'ring')?.ct, p = r && build3D(OCL, rdMol(r.product.molfile)), now = p && stereo(OCL, p.mol, p.given).find(x => x.kind === 'ring')?.ct;
    if (was !== 'cis' || now !== 'trans') bad.push(`ring ${was} → ${now}`);
  }
  if (!bad.length) structuredClone({ hal, r });
  if (bad.length) fail++;
  console.log(bad.length ? 'FAIL' : 'ok  ', 'SN2', name.padEnd(62), rate ?? 'no verdict', r ? rdMol(r.product.molfile) : '', r?.product.rs ?? '', bad.length ? `(${bad.join('; ')})` : '');
}
// No C–Cl, C–Br or C–I on an sp³ carbon: why
for (const [name, smi, want] of [['fluoroethane', 'CCF', 'F'], ['ethanol', 'CCO', 'OH'], ['bromobenzene', 'Brc1ccccc1', 'sp2'], ['vinyl bromide', 'C=CBr', 'sp2'], ['butane', 'CCCC', '']]) {
  const { rows, why } = rxHal(smi).hal, ok = !rows.length && why.join() === want;
  if (!ok) fail++;
  console.log(ok ? 'ok  ' : 'FAIL', 'no C–X', name.padEnd(14), why.join() || '(nothing)', ok ? '' : `(want ${want}, ${rows.length} rows)`);
}
// SN1: [name, SMILES, rate, shift note shown?, flat cation drawn?, relation of the two products ('' = none shown),
// textbook products, index of the C–O carbon in their SMILES]. The products must match the textbook ones (RDKit canonical
// isomeric SMILES), RDKit must agree on the relation (same SMILES, or one is the other with every @/@@ inverted), and each
// product's R/S must be RDKit's. MMFF94s+ builds a ring 3° carbocation at 354.5–354.7° (flat is 360°), under the 355°
// check: it is not drawn, but its two faces still give the two products.
const SN1 = [
  ['(S)-3-bromo-3-methylhexane', 'CCC[C@@](C)(Br)CC', 'good', false, true, 'enantiomers', ['CCC[C@@](C)(O)CC', 'CCC[C@](C)(O)CC'], 3],
  ['2-bromo-2-methylpropane', 'CC(C)(C)Br', 'good', false, true, 'same', ['CC(C)(C)O'], 1],
  ['(3R)-1-bromo-1,3-dimethylcyclohexane', 'C[C@@H]1CCCC(C)(Br)C1', 'good', false, false, 'diastereomers', ['C[C@@H]1CCC[C@](C)(O)C1', 'C[C@@H]1CCC[C@@](C)(O)C1'], 5],
  ['2-bromo-3-methylbutane', 'CC(Br)C(C)C', 'slow', true, true, 'enantiomers', ['C[C@H](O)C(C)C', 'C[C@@H](O)C(C)C'], 1],
  ['2-bromobutane', 'CC(Br)CC', 'slow', false, true, 'enantiomers', ['C[C@H](O)CC', 'C[C@@H](O)CC'], 1],
  ['bromomethane', 'CBr', 'none', false], ['1-bromopropane', 'CCCBr', 'none', false],
  ['benzyl bromide', 'BrCc1ccccc1', 'resonance', false, true, 'same', ['OCc1ccccc1'], 1], ['3-bromoprop-1-ene', 'C=CCBr', 'resonance', false, true, 'same', ['C=CCO'], 2],
  ['2-bromo-3,3-dimethylbutane (methyl shift)', 'CC(Br)C(C)(C)C', 'slow', true, true, 'enantiomers', ['C[C@H](O)C(C)(C)C', 'C[C@@H](O)C(C)(C)C'], 1],
];
if (rdCIP('CCC[C@@](C)(Br)CC').CIP_atoms.map(([, l]) => l).join() !== '(S)') { fail++; console.log('FAIL SN1 test SMILES is not (S)-3-bromo-3-methylhexane'); }
if (rdCIP('C[C@@H]1CCCC(C)(Br)C1').CIP_atoms.find(([i]) => i === 1)?.[1] !== '(R)') { fail++; console.log('FAIL SN1 test SMILES is not (3R)-1-bromo-1,3-dimethylcyclohexane'); }
for (const [name, smi, rate, shift, flat, rel, want = [], at] of SN1) {
  const { mol, hal } = rxHal(smi), h = hal.rows[0], bad = [], t = performance.now();
  if (sn1Rate(h) !== rate) bad.push(`rate ${sn1Rate(h)}`);
  if ((!['none', 'resonance'].includes(rate) && !!h.shift) !== shift) bad.push(`shift ${h.shift}`); // the page shows it only then
  const r = rate === 'none' ? null : sn1(OCL, mol, h.c, h.x);
  if (r) {
    if (!r.cation || !!r.cation.flat !== flat) bad.push(`cation ${r.cation ? `${r.cation.sum.toFixed(1)}°` : 'not built'}`);
    if (!r.products) bad.push('no products');
    else {
      const got = r.products.map(p => rdMol(p.molfile)), exp = want.map(rdCan);
      if (got.some(g => !exp.includes(g)) || exp.some(e => !got.includes(e))) bad.push(`products ${got.join(' ')}`);
      const [a, b] = got, rd = rdCan(a) === rdCan(b) ? 'same' : rdCan(a) === rdCan(b.replace(/@@?/g, x => (x === '@' ? '@@' : '@'))) ? 'enantiomers' : 'diastereomers';
      if (r.rel !== rel || rd !== rel) bad.push(`relation ${r.rel} (RDKit ${rd})`);
      r.products.forEach((p, i) => { const w = want[exp.indexOf(got[i])], l = w && (rdCIP(w).CIP_atoms.find(([j]) => j === at)?.[1][1] ?? null);
        if (w && (p.rs ?? null) !== l) bad.push(`R/S ${p.rs} (RDKit ${l})`); });
      structuredClone(r);
    }
  }
  if (bad.length) fail++;
  console.log(bad.length ? 'FAIL' : 'ok  ', 'SN1', name.padEnd(42), rate, r?.cation ? `${r.cation.sum.toFixed(1)}°` : '', r?.rel ?? '', r?.products?.map(p => p.rs).join(' ') ?? '', `${(performance.now() - t).toFixed(0)}ms`, bad.length ? `(${bad.join('; ')})` : '');
}
// E2: [name, SMILES, textbook alkenes by β carbon (heavy-atom index: { SMILES of what its H give: 'all' | 'one' | 'none' }),
// Zaitsev alkene ('' when all alkenes are alike substituted), can it form?, chair: the X axial or equatorial in the model]. Each H's alkene must match RDKit reading
// the molfile (E/Z from the 3D), and its E/Z must be RDKit's CIP label. Menthyl chloride: only 2-menthene (toward CH₃),
// against Zaitsev; neomenthyl chloride: 2- and 3-menthene, 3-menthene the Zaitsev one.
const MENTHYL = 'C[C@@H]1CC[C@H]([C@@H](C1)Cl)C(C)C', NEOMENTHYL = 'C[C@@H]1CC[C@H]([C@H](C1)Cl)C(C)C';
const M2 = 'C[C@@H]1CC[C@H](C=C1)C(C)C', M3 = 'C[C@@H]1CCC(=CC1)C(C)C'; // 2-menthene (C2=C3), 3-menthene (C3=C4)
const E2 = [
  ['(R)-2-bromobutane', 'C[C@@H](Br)CC', { 0: { 'C=CCC': 'all' }, 3: { 'C/C=C/C': 'one', 'C/C=C\\C': 'one' } }, 'C/C=C/C C/C=C\\C', true],
  ['meso-1,2-dibromo-1,2-diphenylethane', 'Br[C@H](c1ccccc1)[C@@H](Br)c1ccccc1', { 8: { 'Br/C(c1ccccc1)=C/c1ccccc1': 'all' } }, '', true],
  ['(1R,2R)-1,2-dibromo-1,2-diphenylethane', 'Br[C@H](c1ccccc1)[C@H](Br)c1ccccc1', { 8: { 'Br/C(c1ccccc1)=C\\c1ccccc1': 'all' } }, '', true],
  ['menthyl chloride', MENTHYL, { 4: { [M3]: 'none' }, 6: { [M2]: 'one' } }, M3, false, 'equatorial'],
  ['neomenthyl chloride', NEOMENTHYL, { 4: { [M3]: 'all' }, 6: { [M2]: 'one' } }, M3, true, 'axial'],
  ['2-bromo-2-methylbutane', 'CCC(C)(C)Br', { 1: { 'CC=C(C)C': 'all' }, 3: { 'C=C(C)CC': 'all' }, 4: { 'C=C(C)CC': 'all' } }, 'CC=C(C)C', true],
];
for (const s of ['Br[C@H](c1ccccc1)[C@@H](Br)c1ccccc1', 'Br[C@H](c1ccccc1)[C@H](Br)c1ccccc1']) // meso: its own mirror image; (1R,2R)
  if ((rdCan(s) === rdCan(s.replace(/@@?/g, x => (x === '@' ? '@@' : '@')))) !== s.includes('@@') || (!s.includes('@@') && rdCIP(s).CIP_atoms.map(([, l]) => l).join() !== '(R),(R)')) { fail++; console.log('FAIL E2 test SMILES', s); }
for (const [name, smi, want, zaitsev, zOk, xChair] of E2) {
  const { mol, hal, atoms } = rxHal(smi), h = hal.rows[0], ch = chairs(mol, atoms), t = performance.now(), r = e2(OCL, mol, atoms, ch, h.c, h.x), bad = [];
  const molOf = i => r.hs.find(q => q.alkene === i)?.molfile, rdEZ = mf => { const m = RDKit.get_mol(mf, JSON.stringify({ removeHs: true })), b = JSON.parse(m.get_stereo_tags()).CIP_bonds; m.delete(); return b[0]?.[2]?.[1] ?? null; };
  for (const [beta, outs] of Object.entries(want)) {
    const qs = r.hs.filter(q => r.betas[q.beta] === +beta);
    if (!qs.length) { bad.push(`no H on β ${beta}`); continue; }
    for (const [w, how] of Object.entries(outs)) {
      const give = qs.filter(q => q.alkene !== null && rdMol(molOf(q.alkene)) === rdCan(w)), okGive = give.filter(q => q.ok);
      const n = how === 'all' ? qs.length : how === 'one' ? 1 : 0;
      if (okGive.length !== n || (how === 'all' && give.length !== qs.length)) bad.push(`β ${beta} → ${w}: ${okGive.length} of ${qs.length} H`);
    }
  }
  r.alkenes.forEach((a, i) => { const ez = rdEZ(molOf(i)); if (a.ez !== ez) bad.push(`alkene ${rdMol(molOf(i))} E/Z ${a.ez} (RDKit ${ez})`); });
  const z = r.alkenes.map((a, i) => (a.zaitsev ? rdMol(molOf(i)) : null)).filter(Boolean).sort().join(' ');
  if (z !== (zaitsev ? zaitsev.split(' ').map(rdCan).sort().join(' ') : '')) bad.push(`Zaitsev ${z}`); // '' = one kind of alkene only
  if (r.alkenes.filter(a => a.zaitsev).some(a => a.possible !== zOk)) bad.push(`Zaitsev ${zOk ? 'cannot' : 'can'} form`);
  if (xChair) { const s = ch.find(c => c.atoms.includes(h.c))?.subs.find(q => q.sub === h.x);
    if (!s || (s.ax ? 'axial' : 'equatorial') !== xChair) bad.push(`Cl ${s ? (s.ax ? 'axial' : 'equatorial') : 'not on a chair'}`);
    if (xChair === 'equatorial' && ch[0].subs.some(q => !q.h && q.ax)) bad.push('menthyl chloride not all equatorial'); }
  if (!bad.length) structuredClone(r);
  if (bad.length) fail++;
  console.log(bad.length ? 'FAIL' : 'ok  ', 'E2', name.padEnd(40), r.alkenes.map((a, i) => `${rdMol(molOf(i))}${a.ez ? ` ${a.ez}` : ''}${a.possible ? '' : ' (no anti H)'}${a.zaitsev ? ' Zaitsev' : ''}`).join(', '), `${(performance.now() - t).toFixed(0)}ms`, bad.length ? `(${bad.join('; ')})` : '');
}
// Addition to a C=C: [alkene, SMILES, reagent, kind, textbook products]. The products from the two faces must be the
// textbook ones (RDKit reading the molfiles). The textbook SMILES are checked on their own with RDKit: a pair is
// (R,R)/(S,S) (same CIP label at both carbons) and mirror images, a meso compound is its own mirror image with R and S.
const mirrorS = s => s.replace(/@@?/g, x => (x === '@' ? '@@' : '@'));
const RR = s => [s, mirrorS(s)], BR = 'C[C@@H](Br)[C@@H](C)Br', DIOL = 'C[C@@H](O)[C@@H](C)O', TBR = 'Br[C@H]1CCCC[C@@H]1Br';
const MESO_BR = 'C[C@@H](Br)[C@H](C)Br', MESO_DIOL = 'C[C@@H](O)[C@H](C)O', CIS_DIOL = 'O[C@H]1CCCC[C@H]1O';
for (const [s, meso] of [[BR, false], [DIOL, false], [TBR, false], [MESO_BR, true], [MESO_DIOL, true], [CIS_DIOL, true], [SMI['cis-1,2-dimethylcyclohexane'], true]]) {
  const l = rdCIP(s).CIP_atoms.map(([, x]) => x), self = rdCan(s) === rdCan(mirrorS(s));
  if (l.length !== 2 || (l[0] === l[1]) === meso || self !== meso) { fail++; console.log('FAIL addition test SMILES', s, l.join(), self ? 'its own mirror image' : 'chiral'); }
}
const ADD = [
  ['cis-2-butene', SMI['cis-2-butene'], 'br2', 'enantiomers', RR(BR)], ['trans-2-butene', SMI['trans-2-butene'], 'br2', 'meso', [MESO_BR]],
  ['cis-2-butene', SMI['cis-2-butene'], 'oso4', 'meso', [MESO_DIOL]], ['trans-2-butene', SMI['trans-2-butene'], 'oso4', 'enantiomers', RR(DIOL)],
  ['cyclohexene', 'C1CCC=CC1', 'br2', 'enantiomers', RR(TBR)], ['cyclohexene', 'C1CCC=CC1', 'oso4', 'meso', [CIS_DIOL]],
  ['1,2-dimethylcyclohexene', 'CC1=C(C)CCCC1', 'h2', 'meso', [SMI['cis-1,2-dimethylcyclohexane']]],
  ['propene', 'C=CC', 'br2', 'enantiomers', ['C[C@@H](Br)CBr', 'C[C@H](Br)CBr']], ['propene', 'C=CC', 'oso4', 'enantiomers', ['C[C@@H](O)CO', 'C[C@H](O)CO']],
  ['2-methylpropene', 'C=C(C)C', 'br2', 'none', ['CC(C)(Br)CBr']], ['cis-2-butene', SMI['cis-2-butene'], 'h2', 'none', ['CCCC']],
];
for (const [name, smi, rg, kind, want] of ADD) {
  const { mol } = build3D(OCL, smi), d = doubleBonds(OCL, mol), t = performance.now(), r = d.rows.length === 1 && addition(OCL, mol, ...d.rows[0].atoms, rg), bad = [];
  if (!r?.products) bad.push(d.rows.length !== 1 ? `${d.rows.length} C=C rows` : 'no products');
  else {
    const got = [...new Set(r.products.map(p => rdMol(p.molfile)))].sort(), exp = [...new Set(want.map(rdCan))].sort();
    if (got.join(' ') !== exp.join(' ')) bad.push(`products ${got.join(' ')} (want ${exp.join(' ')})`);
    if (r.kind !== kind) bad.push(`kind ${r.kind}`);
    structuredClone(r);
  }
  if (bad.length) fail++;
  console.log(bad.length ? 'FAIL' : 'ok  ', 'add', `${name} + ${rg}`.padEnd(32), r?.kind ?? '', r?.products ? [...new Set(r.products.map(p => rdMol(p.molfile)))].join(' ') : '', `${(performance.now() - t).toFixed(0)}ms`, bad.length ? `(${bad.join('; ')})` : '');
}
for (const [smi, rows, arom] of [['c1ccccc1', 0, true], ['C=Cc1ccccc1', 1, true], ['C=C=C', 0, false], ['C=CC=C', 1, false], ['C=CCC=C', 1, false], ['CC=CCC=C', 2, false]]) {
  const d = doubleBonds(OCL, build3D(OCL, smi).mol);
  if (d.rows.length !== rows || d.aromatic !== arom) { fail++; console.log('FAIL C=C rows', smi, d.rows.length, d.aromatic, `(want ${rows}, ${arom})`); }
}
// Equivalent C–X bonds are listed once (meso-2,3-dibromobutane: mirror images; 1,2-dibromoethane), different ones apart
for (const [smi, want] of [['C[C@@H](Br)[C@H](C)Br', 1], ['BrCCBr', 1], ['BrCCCl', 2], ['C[C@@H](Br)[C@@H](C)Br', 1]]) {
  const n = rxHal(smi).hal.rows.length;
  if (n !== want) { fail++; console.log('FAIL C–X rows', smi, n, `(want ${want})`); }
}

// worker.js turns flipped-chair coordinates into the molfile's frame, which writes y and z negated: check that rule.
{ const { mol } = build3D(OCL, 'CC1CCCCC1'), lines = mol.toMolfile().split('\n');
  const bad = Array.from({ length: mol.getAllAtoms() }, (_, i) => lines[4 + i].slice(0, 30).trim().split(/\s+/).map(Number))
    .some(([x, y, z], i) => Math.abs(x - mol.getAtomX(i)) > 1e-3 || Math.abs(y + mol.getAtomY(i)) > 1e-3 || Math.abs(z + mol.getAtomZ(i)) > 1e-3);
  if (bad) { fail++; console.log('FAIL molfile frame is no longer (x, −y, −z): fix the turn in worker.js'); } }

// The worker posts these to the page, so they must survive structured cloning (OCL's own arrays do not).
try { const m = build3D(OCL, 'CC1CCC2CCCCC2C1').mol, at = atomInfo(m); structuredClone({ b: newmanBonds(OCL, m, at), r: ringAngles(m, at), c: chairs(m, at) });
  const n = build3D(OCL, 'c1ccc2[nH]ccc2c1').mol; structuredClone({ p: polarity(n), a: aromaticRings(n), e: equivalentH(n) }); }
catch (e) { fail++; console.log('FAIL Shape data cannot be posted from the worker:', e.message); }

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
