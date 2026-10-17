// Self-check: `npm install && node test.mjs`. Builds each molecule in 3D and checks the functional groups found.
import OCL from 'openchemlib';
import initRDKit from '@rdkit/rdkit';
import { build3D, findGroups, atomInfo, kitParts, stereo, cipOrder, superimposable, newmanBonds, conformer, ringAngles, stericPairs, chairs, diaxial, flipChair } from './chem.js';
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

// The worker posts these to the page, so they must survive structured cloning (OCL's own arrays do not).
try { const m = build3D(OCL, 'CC1CCC2CCCCC2C1').mol, at = atomInfo(m); structuredClone({ b: newmanBonds(OCL, m, at), r: ringAngles(m, at), c: chairs(m, at) }); }
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
