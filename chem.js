// Functional groups. Each SMARTS lists the atoms to highlight first;
// `n` = how many leading (non-excluded) pattern atoms get highlighted.
// Written for OpenChemLib's SMARTS parser: negative recursion ($) only.
export const GROUPS = [
  { id: 'alkene', name: 'Alkene', smarts: ['[CX3;!a]=[CX3;!a]'], n: 2,
    info: 'A C=C double bond. Both carbons are sp² and flat, and the double bond cannot rotate, so cis/trans isomers are possible.' },
  { id: 'alkyne', name: 'Alkyne', smarts: ['[CX2]#[CX2]'], n: 2,
    info: 'A C≡C triple bond. Both carbons are sp, so the atoms around them sit on a straight line (180°).' },
  { id: 'aromatic', name: 'Aromatic ring', smarts: [], rings: true, // from OCL's ring perception: its SMARTS misses fused 5-rings (indole)
    info: 'A flat ring with 4n+2 π electrons shared around the whole ring, like benzene. The ring is stable, so it reacts by substitution, not addition.' },
  { id: 'alkylhalide', name: 'Alkyl halide', smarts: ['[CX4][F,Cl,Br,I]'], n: 2,
    info: 'A halogen (F, Cl, Br, I) on an sp³ carbon. The C–X bond is polar, so the carbon is a target for nucleophiles (SN1, SN2) and for elimination (E1, E2).' },
  { id: 'arylhalide', name: 'Aryl halide', smarts: ['c[F,Cl,Br,I]'], n: 2,
    info: 'A halogen bonded straight to an aromatic ring carbon. Unlike an alkyl halide, it does not react by SN1 or SN2.' },
  { id: 'alcohol', name: 'Alcohol', smarts: ['[OX2H][CX4]'], n: 2,
    info: 'An –OH on an sp³ carbon. It gives and accepts hydrogen bonds, so small alcohols mix with water and boil higher than alkanes of the same size.' },
  { id: 'phenol', name: 'Phenol', smarts: ['[OX2H]c'], n: 2,
    info: 'An –OH bonded straight to an aromatic ring. It is more acidic than an alcohol (pKa ≈ 10, not ≈ 16) because the ring spreads out the negative charge of the phenoxide ion.' },
  { id: 'ether', name: 'Ether', smarts: ['[OX2;!a]([#6;!$(C=[O,S,N])])[#6;!$(C=[O,S,N])]'], n: 3,
    info: 'An oxygen between two carbons (C–O–C). It accepts hydrogen bonds but has no O–H to give one. Ethers are fairly unreactive and make common solvents.' },
  { id: 'thiol', name: 'Thiol', smarts: ['[SX2H][#6]'], n: 2,
    info: 'An –SH group, the sulfur version of an alcohol. It is more acidic than an alcohol, and two thiols oxidize to a disulfide (S–S), as in cysteine.' },
  { id: 'sulfide', name: 'Sulfide', smarts: ['[SX2;!a]([#6;!$(C=[O,S,N])])[#6;!$(C=[O,S,N])]'], n: 3,
    info: 'A sulfur between two carbons (C–S–C), the sulfur version of an ether. Methionine has one.' },
  { id: 'aldehyde', name: 'Aldehyde', smarts: ['[CX3H1](=O)[#6]', '[CX3H2]=O'], n: 2, carbonH: true,
    info: 'A carbonyl (C=O) with at least one H on the carbonyl carbon. Nucleophiles attack the C=O carbon, and it oxidizes easily to a carboxylic acid.' },
  { id: 'ketone', name: 'Ketone', smarts: ['[CX3](=O)([#6])[#6]'], n: 2,
    info: 'A carbonyl (C=O) with two carbons on the carbonyl carbon. Nucleophiles attack it like an aldehyde, but it resists oxidation.' },
  { id: 'acid', name: 'Carboxylic acid', smarts: ['[CX3](=O)[OX2H]'], n: 3,
    info: 'A –COOH group: a C=O and an –OH on the same carbon. It is a weak acid (pKa ≈ 4–5) and loses H⁺ to form a carboxylate ion.' },
  { id: 'ester', name: 'Ester', smarts: ['[CX3](=O)([OX2][#6;!$(C=[O,S,N])])[#6]', '[CX3H1](=O)[OX2][#6;!$(C=[O,S,N])]'], n: 3,
    info: 'A –COO– link between two carbon groups. A carboxylic acid and an alcohol make it, and hydrolysis splits it back. Many esters smell fruity.' },
  { id: 'amide', name: 'Amide', smarts: ['[#6X3](=O)[#7X3]'], n: 3,
    info: 'A carbonyl bonded to nitrogen (C(=O)–N). Resonance gives the C–N bond partial double-bond character, so the group is flat and the N is not basic. Peptide bonds are amides.' },
  { id: 'acylhalide', name: 'Acyl halide', smarts: ['[CX3](=O)[F,Cl,Br,I]'], n: 3,
    info: 'A carbonyl bonded to a halogen (–COCl). It is the most reactive carboxylic acid derivative and reacts fast with water, alcohols and amines.' },
  { id: 'anhydride', name: 'Acid anhydride', smarts: ['[CX3](=O)[OX2][CX3]=O'], n: 5,
    info: 'Two acyl groups joined through one oxygen (–CO–O–CO–). It is a strong acylating agent; acetic anhydride makes aspirin from salicylic acid.' },
  { id: 'amine', name: 'Amine', smarts: ['[NX3;!a;!$(N~[!#6;!#1]);!$(N[#6]=[O,S,N]);!$(N=*)]'], n: 1,
    info: 'A nitrogen bonded only to carbons or hydrogens, with no C=O next to it. Its lone pair makes it a base and a nucleophile.' },
  { id: 'imine', name: 'Imine', smarts: ['[CX3;!a;!$(C-[!#6;!#1])]=[NX2;!a;!$(N-[O,N])]'], n: 2,
    info: 'A C=N double bond. An aldehyde or ketone forms one when it reacts with a primary amine.' },
  { id: 'nitrile', name: 'Nitrile', smarts: ['[CX2]#[NX1]'], n: 2,
    info: 'A –C≡N group. The carbon is sp and linear. Hydrolysis turns a nitrile into a carboxylic acid.' },
  { id: 'nitro', name: 'Nitro', smarts: ['[NX3+](=O)([O-])[#6]', '[NX3](=O)(=O)[#6]'], n: 3,
    info: 'An –NO₂ group on carbon. It pulls electrons away strongly; reducing it gives an amine (nitrobenzene → aniline).' },
];

// SMILES → OCL Molecule with explicit H and 3D coordinates (OCL torsion-library conformer, then MMFF94s+).
// Salts: only the largest fragment is kept.
// `given` = stereo parities from the SMILES itself (3 = left open), keyed by heavy-atom index, which the
// 3D molecule keeps (it only appends H). The 3D model always picks one isomer, so this tells which picks were set.
export function build3D(OCL, smiles) {
  // getFragments, not stripSmallFragments: the latter neutralizes charges (acetate would become acetic acid)
  const frags = OCL.Molecule.fromSmiles(smiles).getFragments();
  const mol = frags.reduce((a, b) => (b.getAtoms() > a.getAtoms() ? b : a));
  const stripped = frags.length > 1;
  mol.ensureHelperArrays(OCL.Molecule.cHelperParities);
  const given = { atoms: [], bonds: {} };
  for (let a = 0; a < mol.getAtoms(); a++) given.atoms.push(mol.getAtomParity(a));
  for (let b = 0; b < mol.getBonds(); b++) given.bonds[pairKey(mol.getBondAtom(0, b), mol.getBondAtom(1, b))] = mol.getBondParity(b);
  const m3 = new OCL.ConformerGenerator(42).getOneConformerAsMolecule(mol);
  if (!m3) throw new Error('Could not build a 3D model for this structure.');
  let minimized = true;
  try { new OCL.ForceFieldMMFF94(m3, 'MMFF94s+').minimise(); }
  catch { minimized = false; } // atom types MMFF does not cover (e.g. metals): keep the raw conformer
  return { mol: m3, stripped, minimized, given };
}

const pairKey = (a, b) => (a < b ? `${a}-${b}` : `${b}-${a}`);
const nbrs = (mol, a) => Array.from({ length: mol.getAllConnAtoms(a) }, (_, k) => mol.getConnAtom(a, k));
const pos = (mol, a) => [mol.getAtomX(a), mol.getAtomY(a), mol.getAtomZ(a)];
const sub = (p, q) => [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
const dot = (p, q) => p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
const cross = (p, q) => [p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]];
const angle = (p, q) => Math.acos(Math.max(-1, Math.min(1, dot(p, q) / Math.sqrt(dot(p, p) * dot(q, q))))) * 180 / Math.PI;

// ---- Structure tab: hybridization, lone pairs, bond angles ----
// Bond orders below are OCL's Kekulé orders (aromatic rings come back as alternating 1/2).
const VALENCE_E = { H: 1, B: 3, C: 4, N: 5, O: 6, F: 7, Si: 4, P: 5, S: 6, Cl: 7, Se: 6, Br: 7, I: 7 };
export const IDEAL = { sp3: 109.5, sp2: 120, sp: 180 };
const SECOND_ROW = new Set([6, 7, 8]); // C, N, O: textbooks assign hybridization to these

function piBonds(mol, a) {
  let n = 0;
  for (let k = 0; k < mol.getAllConnAtoms(a); k++) n += mol.getConnBondOrder(a, k) - 1;
  return n;
}

// (valence electrons − charge − electrons used in bonds) / 2; null when that is not a whole number (radicals, metals).
export function lonePairs(mol, a) {
  const v = VALENCE_E[mol.getAtomLabel(a)];
  if (v === undefined) return null;
  let used = 0;
  for (let k = 0; k < mol.getAllConnAtoms(a); k++) used += mol.getConnBondOrder(a, k);
  const lp = (v - mol.getAtomCharge(a) - used) / 2;
  return lp >= 0 && Number.isInteger(lp) ? lp : null;
}

// Textbook rule: count π bonds first (2 → sp, 1 → sp2); otherwise count σ bonds + lone pairs (4 → sp3, 3 → sp2).
// An sp3-by-count C, N or O whose lone pair sits next to a π bond (or an empty p orbital) shares it by resonance,
// so it counts as sp2 (amide N, aniline N, pyrrole N, ester and phenol O).
// shortcut: geometry that blocks the overlap (bridgehead amides) is not checked; add a ring test if one comes up.
export function hybridization(mol, a) {
  if (!SECOND_ROW.has(mol.getAtomicNo(a))) return null;
  const lp = lonePairs(mol, a);
  if (lp === null) return null;
  const pi = piBonds(mol, a);
  if (pi >= 2) return { hyb: 'sp' };
  if (pi === 1) return { hyb: 'sp2' };
  const groups = mol.getAllConnAtoms(a) + lp;
  if (groups === 3) return { hyb: 'sp2' }; // e.g. a carbocation: 3 σ bonds and an empty p orbital
  if (groups !== 4) return null;
  const conjugated = lp > 0 && nbrs(mol, a).some(b => SECOND_ROW.has(mol.getAtomicNo(b))
    && (piBonds(mol, b) > 0 || (mol.getAllConnAtoms(b) === 3 && lonePairs(mol, b) === 0)));
  return conjugated ? { hyb: 'sp2', resonance: true } : { hyb: 'sp3' };
}

// Per atom: element, hybridization, lone pairs, neighbors, and every bond angle at the atom (from the 3D model).
export function atomInfo(mol) {
  const out = [];
  for (let a = 0; a < mol.getAllAtoms(); a++) {
    const nb = nbrs(mol, a), p = pos(mol, a), angles = [];
    for (let i = 0; i < nb.length; i++) for (let j = i + 1; j < nb.length; j++)
      angles.push([nb[i], nb[j], angle(sub(pos(mol, nb[i]), p), sub(pos(mol, nb[j]), p))]);
    out.push({ el: mol.getAtomLabel(a), charge: mol.getAtomCharge(a), lp: lonePairs(mol, a), nbrs: nb, angles, ...hybridization(mol, a) });
  }
  return out;
}

// ---- Kit tab: Darling Molecular Visions parts (rules from Dan Berger's Bluffton University guide) ----
// Aromatic rings are counted in one Kekulé structure. Anything the guide does not cover goes to `check`.
export const KIT_PARTS = [
  { id: 'black', label: 'Black tetrahedral pieces (C)' },
  { id: 'blue', label: 'Blue tetrahedral pieces (N)' },
  { id: 'red', label: 'Red tetrahedral pieces (O)' },
  { id: 'dbl', label: 'Gray double-bond pieces (C=C, C=N)' },
  { id: 'half', label: 'Gray + red half-bond pairs (C=O)' },
  { id: 'end', label: 'Gray end-pieces (one per double-bonded C)' },
  { id: 'tri', label: 'Gray triple-bond pieces (C≡C)' },
  { id: 'bend', label: 'Ring atoms to build with bendable pieces (3- and 4-membered rings)' },
  { id: 'h', label: 'White balls for H (optional: bare arm ends work too)' },
  { id: 'lp', label: 'Colored balls for lone pairs (optional)' },
];

const KIT_BONDS = { 'C=C': 'dbl', 'C=N': 'dbl', 'C=O': 'half', 'C≡C': 'tri' };

export function kitParts(mol) {
  const n = Object.fromEntries(KIT_PARTS.map(p => [p.id, 0]));
  n.redOne = 0; // red pieces if each 2-bond O is one piece instead of two (the guide allows both)
  const check = {};
  const flag = what => { check[what] = (check[what] || 0) + 1; };
  const el = a => mol.getAtomLabel(a);
  const small = a => mol.isRingAtom(a) && mol.getAtomRingSize(a) <= 4;
  const smallBond = b => mol.isRingBond(b) && mol.getBondRingSize(b) <= 4;
  const bondName = b => {
    const [x, y] = [el(mol.getBondAtom(0, b)), el(mol.getBondAtom(1, b))].sort();
    return `${x}${mol.getBondOrder(b) === 2 ? '=' : '≡'}${y}`;
  };

  for (let b = 0; b < mol.getAllBonds(); b++) {
    if (mol.getBondOrder(b) < 2) continue;
    const name = bondName(b);
    if (smallBond(b)) flag('Double or triple bond in a 3- or 4-membered ring: check your kit');
    else if (KIT_BONDS[name]) n[KIT_BONDS[name]]++;
    else if (name === 'C≡N') flag('C≡N (nitrile): this kit cannot show it');
    else flag(`${name} bond: check your kit`);
  }

  for (let a = 0; a < mol.getAllAtoms(); a++) {
    const e = el(a);
    if (e === 'H') { n.h++; continue; }
    if (!['C', 'N', 'O'].includes(e)) { flag(`${e} atoms: check your kit`); continue; }
    const conn = mol.getAllConnAtoms(a), charge = mol.getAtomCharge(a), lp = lonePairs(mol, a);
    const multi = [];
    for (let k = 0; k < conn; k++) if (mol.getConnBondOrder(a, k) > 1) multi.push(mol.getConnBond(a, k));
    if (!multi.length) {
      if (conn + lp !== 4) flag(`${e} with ${conn} bonds and no lone pair (flat): check your kit`);
      else if (small(a)) n.bend++;
      else {
        n[{ C: 'black', N: 'blue', O: 'red' }[e]] += 2;
        if (e === 'O') n.redOne += conn === 2 ? 1 : 2;
        n.lp += lp;
      }
      continue;
    }
    if (multi.some(b => smallBond(b) || !KIT_BONDS[bondName(b)])) continue; // already flagged with its bond
    if (small(a)) { flag('Double-bonded atom in a 3- or 4-membered ring: check your kit'); continue; }
    const kind = multi.length === 1 ? bondName(multi[0]) : '';
    const ok = !charge && (
      (e === 'C' && kind !== 'C≡C' && kind && conn === 3) ||
      (e === 'C' && kind === 'C≡C' && conn === 2) ||
      (e === 'N' && kind === 'C=N' && conn === 2) ||
      (e === 'O' && kind === 'C=O' && conn === 1));
    if (!ok) flag(`${e} with ${multi.length > 1 ? 'two multiple bonds' : charge ? 'a charge next to a multiple bond' : 'an unusual set of bonds'}: check your kit`);
    else if (e === 'C' && kind !== 'C≡C') n.end++;
    else if (e === 'N') n.blue += 1; // the guide: one blue tetrahedral piece on the N side of C=N
  }
  return { n, check: Object.entries(check).map(([what, count]) => ({ what, count })) };
}

// ---- Stereo tab: double-bond E/Z and cis/trans, ring cis/trans ----
// E/Z comes from OCL's CIP rules on the 3D model. cis/trans is read from the 3D geometry.
// `given` (from build3D) says whether the name or SMILES set each one, or the 3D build picked it.
export function stereo(OCL, mol, given) {
  mol.ensureHelperArrays(OCL.Molecule.cHelperCIP);
  const isH = a => mol.getAtomicNo(a) === 1;
  const out = [];

  for (let b = 0; b < mol.getAllBonds(); b++) {
    if (mol.getBondOrder(b) !== 2 || ![1, 2].includes(mol.getBondParity(b))) continue;
    const [p, q] = [mol.getBondAtom(0, b), mol.getBondAtom(1, b)];
    const cip = mol.getBondCIPParity(b);
    const row = { kind: 'bond', atoms: [p, q, ...nbrs(mol, p), ...nbrs(mol, q)].filter((x, i, s) => s.indexOf(x) === i),
      label: `${mol.getAtomLabel(p)}=${mol.getAtomLabel(q)}`, ez: cip === 1 ? 'E' : cip === 2 ? 'Z' : null,
      given: [1, 2].includes(given.bonds[pairKey(p, q)]) };
    // cis/trans only for 1,2-disubstituted C=C: each carbon has one H and one other group
    const side = (c, d) => { const o = nbrs(mol, c).filter(x => x !== d); return o.length === 2 && o.filter(isH).length === 1 ? o.find(x => !isH(x)) : -1; };
    const sp = side(p, q), sq = side(q, p);
    if (mol.getAtomLabel(p) === 'C' && mol.getAtomLabel(q) === 'C' && sp >= 0 && sq >= 0)
      row.ct = dihedral(mol, sp, p, q, sq) < 90 ? 'cis' : 'trans';
    out.push(row);
  }

  // Ring cis/trans: ring atoms carrying one H and one other group; compare which face each group sits on.
  // shortcut: OCL's ring set stops at 7-membered rings; extend if larger rings come up in class.
  const seen = new Set(), rs = mol.getRingSet();
  for (let r = 0; r < rs.getSize(); r++) {
    const ring = rs.getRingAtoms(r), inRing = new Set(ring);
    const c = ring.map(a => pos(mol, a)).reduce((s, p) => s.map((v, i) => v + p[i] / ring.length), [0, 0, 0]);
    let normal = [0, 0, 0]; // Newell-style: sum of cross products around the ring (atoms come in ring order)
    ring.forEach((a, i) => { normal = normal.map((v, k) => v + cross(sub(pos(mol, a), c), sub(pos(mol, ring[(i + 1) % ring.length]), c))[k]); });
    const subs = [];
    for (const a of ring) {
      const exo = nbrs(mol, a).filter(x => !inRing.has(x));
      if (exo.length !== 2 || exo.filter(isH).length !== 1) continue;
      const g = exo.find(x => !isH(x)), h = exo.find(isH);
      subs.push({ a, g, h, up: dot(normal, sub(pos(mol, g), pos(mol, h))) > 0 });
    }
    for (let i = 0; i < subs.length; i++) for (let j = i + 1; j < subs.length; j++) {
      const [s, t] = [subs[i], subs[j]], key = pairKey(s.a, t.a);
      if (seen.has(key) || !mol.getAtomParity(s.a) || !mol.getAtomParity(t.a)) continue; // not a stereo unit
      seen.add(key);
      out.push({ kind: 'ring', size: ring.length, atoms: [s.a, s.g, s.h, t.a, t.g, t.h], groups: [s.g, t.g],
        ct: s.up === t.up ? 'cis' : 'trans', given: [1, 2].includes(given.atoms[s.a]) && [1, 2].includes(given.atoms[t.a]) });
    }
  }
  return out;
}

function dihedral(mol, a, b, c, d) {
  const [pa, pb, pc, pd] = [a, b, c, d].map(x => pos(mol, x));
  const n1 = cross(sub(pb, pa), sub(pc, pb)), n2 = cross(sub(pc, pb), sub(pd, pc));
  return angle(n1, n2);
}

// Returns [{group, matches: [[atomIndex...], ...], atoms}] for groups found in `mol` (an OCL Molecule).
// `atoms` = all highlighted atoms plus their H (H on heteroatoms; on carbon only for aldehydes).
export function findGroups(OCL, mol) {
  const searcher = new OCL.SSSearcher();
  const out = [];
  for (const g of GROUPS) {
    const seen = new Set(), matches = [];
    if (g.rings) {
      const rs = mol.getRingSet();
      for (let r = 0; r < rs.getSize(); r++) if (rs.isAromatic(r)) matches.push(rs.getRingAtoms(r));
    }
    for (const s of g.smarts) {
      const frag = new OCL.SmilesParser({ smartsMode: 'smarts' }).parseMolecule(s);
      frag.setFragment(true);
      searcher.setMol(frag, mol);
      if (!searcher.findFragmentInMolecule({ countMode: 'overlapping' })) continue;
      for (const m of searcher.getMatchList()) {
        const atoms = m.filter(i => i >= 0).slice(0, g.n);
        const key = [...atoms].sort((a, b) => a - b).join();
        if (!seen.has(key)) { seen.add(key); matches.push(atoms); }
      }
    }
    if (!matches.length) continue;
    const atoms = new Set(matches.flat());
    for (const a of [...atoms]) if (g.carbonH || mol.getAtomicNo(a) !== 6)
      for (let k = 0; k < mol.getAllConnAtoms(a); k++) {
        const b = mol.getConnAtom(a, k);
        if (mol.getAtomicNo(b) === 1) atoms.add(b);
      }
    out.push({ group: g, matches, atoms: [...atoms] });
  }
  return out;
}
