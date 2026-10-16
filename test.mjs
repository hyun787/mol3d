// Self-check: `npm install && node test.mjs`. Builds each molecule in 3D and checks the functional groups found.
import OCL from 'openchemlib';
import { build3D, findGroups } from './chem.js';
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
