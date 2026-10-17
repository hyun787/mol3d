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
  // The conformer generator sometimes ignores a set parity (cis-1,2-dichloroethene comes out trans with seed 42), so a
  // conformer whose 3D parities differ from the ones the SMILES set is thrown away and the next seed tried.
  const keeps = m => {
    const c = m.getCompactCopy();
    c.ensureHelperArrays(OCL.Molecule.cHelperParities);
    if (given.atoms.some((p, a) => [1, 2].includes(p) && c.getAtomParity(a) !== p)) return false;
    for (let b = 0; b < c.getBonds(); b++) {
      const p = given.bonds[pairKey(c.getBondAtom(0, b), c.getBondAtom(1, b))];
      if ([1, 2].includes(p) && c.getBondParity(b) !== p) return false;
    }
    return true;
  };
  let m3 = null;
  for (const seed of [42, ...Array.from({ length: 20 }, (_, i) => i + 1)]) {
    const m = new OCL.ConformerGenerator(seed).getOneConformerAsMolecule(mol);
    if (m && (m3 ??= m) && keeps(m)) { m3 = m; break; }
  }
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

// ---- Stereo tab: R/S stereocenters, double-bond E/Z and cis/trans, ring cis/trans ----
// R/S and E/Z come from OCL's CIP rules on the 3D model. cis/trans is read from the 3D geometry.
// `given` (from build3D) says whether the name or SMILES set each one, or the 3D build picked it.
export function stereo(OCL, mol, given) {
  mol.ensureHelperArrays(OCL.Molecule.cHelperCIP);
  const isH = a => mol.getAtomicNo(a) === 1;
  const out = [];

  // Stereocenters with 4 neighbors (3-neighbor N, S, P and the allene middle carbon are left out).
  // A center is listed only when our own CIP ranking (cipOrder) separates all 4 groups AND gives the same R/S as OCL
  // on this 3D model. That drops pseudoasymmetric centers (r/s: two groups differ only in configuration, which
  // OCL sometimes labels R/S anyway) and the rare centers that need CIP rules 3–5.
  for (let a = 0; a < mol.getAtoms(); a++) {
    const cip = mol.getAtomCIPParity(a), nb = nbrs(mol, a);
    if (![1, 2].includes(cip) || mol.isAtomParityPseudo(a) || nb.length !== 4) continue;
    const order = cipOrder(mol, a), rs = cip === 1 ? 'R' : 'S';
    if (!order || handedness(mol, a, order) !== rs) continue;
    out.push({ kind: 'center', atom: a, atoms: [a, ...nb], rs, order, sets: order.map(b => firstSphere(mol, b, a)), given: [1, 2].includes(given.atoms[a]) });
  }

  for (let b = 0; b < mol.getAllBonds(); b++) {
    if (mol.getBondOrder(b) !== 2 || ![1, 2].includes(mol.getBondParity(b))) continue;
    if (mol.isRingBond(b) && mol.getBondRingSize(b) < 8) continue; // OCL gives ring C=C a parity, but below 8 atoms only cis fits
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
    const ring = Array.from(rs.getRingAtoms(r)), inRing = new Set(ring), normal = ringFrame(mol, ring).n;
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

// ---- CIP priorities: rule 1a (atomic number), then rule 2 (mass), on the hierarchical digraph ----
// Each group is explored sphere by sphere from the center, set by set in order of rank. A double bond gives each end
// a duplicate of the other end; a ring closure gives a duplicate of the atom the path returns to. Duplicates carry
// only phantom atoms (atomic number 0), the same as a missing entry here, and so do lone pairs.
// shortcut: rules 1b and 3–5 are not implemented; a tie they would break leaves the center out (cipOrder → null).
// shortcut: aromatic rings use OCL's Kekulé bonds; IUPAC averages duplicates over Kekulé structures, which only
// changes heteroaromatic rings. Add that if a pyridine-type case ever disagrees with RDKit in test.mjs.
const MASS = { 1: 1.008, 6: 12.011, 7: 14.007, 8: 15.999, 9: 18.998, 14: 28.085, 15: 30.974, 16: 32.06, 17: 35.45, 35: 79.904, 53: 126.904 };
const PHANTOM = { phantom: true, kids: [[], [], []] };

// Returns the 4 neighbor atoms of `center` from priority 1 to 4, or null when two of them tie.
export function cipOrder(mol, center) {
  let budget = 100000; // digraph nodes: a huge ring system gives up instead of freezing the page
  const iso = Array.from({ length: mol.getAllAtoms() }, (_, a) => mol.getAtomMass(a)).some(Boolean);
  const z = n => (n.phantom ? 0 : mol.getAtomicNo(n.atom));
  const mass = n => (n.phantom ? 0 : mol.getAtomMass(n.atom) || MASS[z(n)] || z(n));
  const keys = [null, z, n => z(n) + mass(n) / 1000]; // pass 2: same atomic numbers, then mass
  const node = (atom, parent, dup) => {
    if (--budget < 0) throw new Error('digraph too large');
    return { atom, parent, dup, kids: [] };
  };
  const onPath = (n, b) => { for (let p = n; p; p = p.parent) if (p.atom === b) return true; return false; };

  // Children of a node, sorted from highest to lowest rank under rule `r` (pass 2 sorts by rule 1 first).
  const kids = (n, r) => {
    if (n.kids[r]) return n.kids[r];
    let out = n.kids[1] || n.kids[2];
    if (!out) {
      out = [];
      if (!n.dup) for (let k = 0; k < mol.getAllConnAtoms(n.atom); k++) {
        const b = mol.getConnAtom(n.atom, k), order = mol.getConnBondOrder(n.atom, k);
        const back = b === n.parent.atom;
        if (!back) out.push(node(b, n, onPath(n.parent, b)));
        for (let i = 1; i < order; i++) out.push(node(b, n, true));
      }
    }
    return n.kids[r] = [...out].sort((x, y) => rank(y, x, 1) || (r === 2 ? rank(y, x, 2) : 0));
  };

  // > 0 when node a outranks node b. Compares sphere by sphere; within a sphere, set by set in branch order.
  const rank = (a, b, r) => {
    const key = keys[r];
    let d = key(a) - key(b);
    if (d) return d;
    let A = [a], B = [b];
    while (A.length) {
      const nA = [], nB = [];
      for (let i = 0; i < A.length; i++) {
        const ka = kids(A[i], r), kb = kids(B[i], r);
        for (let j = 0; j < Math.max(ka.length, kb.length); j++) {
          const x = ka[j] || PHANTOM, y = kb[j] || PHANTOM;
          if ((d = key(x) - key(y))) return d;
          nA.push(x); nB.push(y);
        }
      }
      A = nA; B = nB;
    }
    return 0;
  };

  try {
    const root = { atom: center, parent: null };
    const ligs = nbrs(mol, center).map(b => node(b, root, false));
    const cmp = (x, y) => rank(y, x, 1) || (iso ? rank(y, x, 2) : 0);
    ligs.sort(cmp);
    return ligs.some((l, i) => i && !cmp(ligs[i - 1], l)) ? null : ligs.map(l => l.atom);
  } catch { return null; }
}

// The atoms on a group's first atom, as CIP compares them: highest atomic number first, a double bond counted twice.
export function firstSphere(mol, b, center) {
  const out = [];
  for (let k = 0; k < mol.getAllConnAtoms(b); k++) {
    const x = mol.getConnAtom(b, k);
    if (x !== center) for (let i = 0; i < mol.getConnBondOrder(b, k); i++) out.push(x);
  }
  return out.sort((x, y) => mol.getAtomicNo(y) - mol.getAtomicNo(x)).map(x => mol.getAtomLabel(x));
}

// R or S from the 3D model: looking with group 4 pointing away, 1 → 2 → 3 clockwise is R.
export function handedness(mol, center, order) {
  const [v1, v2, v3, v4] = order.map(b => pos(mol, b));
  return dot(sub(v1, v4), cross(sub(v2, v4), sub(v3, v4))) < 0 ? 'R' : 'S';
}

// Does the molecule stack on its mirror image? Compares OCL's canonical code (it includes stereo) of the model and of
// its reflection. That judges the configuration, not this one 3D shape, so a flexible achiral molecule counts as
// the same even when this conformer alone would not overlay its mirror image.
export function superimposable(mol) {
  const m = mol.getCompactCopy();
  for (let a = 0; a < m.getAllAtoms(); a++) m.setAtomX(a, -m.getAtomX(a));
  return m.getIDCode() === mol.getIDCode();
}

// ---- Shape tab: Newman projections down C–C bonds, and ring angles ----
const unit = v => { const n = Math.hypot(...v); return v.map(x => x / n); };

// C–C single bonds outside rings with both carbons sp³, one per set of symmetry-equivalent bonds.
// Newman view down the bond: front carbon F (the lower atom index) toward the viewer, back carbon B behind it.
// Angles go around the F→B axis (right-handed), 0° = the front `x` group. `x` = the pair whose dihedral the slider sets:
// the reference groups when each carbon has exactly one non-H group (then 180° is anti), otherwise the first non-H
// (or first H) group on each carbon. `phi0` = that dihedral in the 3D model; `side` = atoms that turn with B.
// `levels` order the six forms at 0°, 60° … 300° (0 = most stable) by the textbook rule: eclipsed is above staggered;
// among eclipsed forms, more non-H/non-H eclipsed pairs is higher; among staggered forms, more gauche pairs of sp³
// carbons is higher. For alkyl groups this is the order the textbook increments give (H/H 4.0, H/CH3 6.0, CH3/CH3
// 11 kJ/mol eclipsed, gauche CH3/CH3 3.8), checked in test.mjs.
// Groups covered: H, sp³ carbon, and a "flat" carbon (sp² with only C and H on it: a C=C or benzene-ring carbon).
// A flat group counts when eclipsed but not when gauche: it turns aside, and MMFF94s+ puts gauche only 0.3–0.9 kJ/mol
// above anti for propylbenzene and 1-pentene (3.3 for butane), also checked in test.mjs. `flat` = the bond has one.
// shortcut: all alkyls count alike, and polar groups are not ranked (1,2-difluoroethane prefers gauche; C=O carbons
// carry dipoles); add them only with checked data.
export function newmanBonds(OCL, mol, atoms) {
  mol.ensureHelperArrays(OCL.Molecule.cHelperSymmetrySimple);
  const isH = a => mol.getAtomicNo(a) === 1;
  const sp3C = a => mol.getAtomicNo(a) === 6 && atoms[a].hyb === 'sp3';
  const flatC = a => mol.getAtomicNo(a) === 6 && atoms[a].hyb === 'sp2' && nbrs(mol, a).every(n => [1, 6].includes(mol.getAtomicNo(n)));
  const seen = new Set(), out = [];
  for (let b = 0; b < mol.getAllBonds(); b++) {
    const [F, B] = [mol.getBondAtom(0, b), mol.getBondAtom(1, b)].sort((p, q) => p - q);
    if (mol.getBondOrder(b) !== 1 || mol.isRingBond(b) || !sp3C(F) || !sp3C(B)) continue;
    const key = [mol.getSymmetryRank(F), mol.getSymmetryRank(B)].sort((p, q) => p - q).join();
    if (seen.has(key)) continue;
    seen.add(key);
    const front = nbrs(mol, F).filter(a => a !== B), back = nbrs(mol, B).filter(a => a !== F);
    const heavy = s => s.filter(a => !isH(a));
    const ref = heavy(front).length === 1 && heavy(back).length === 1;
    const x = [heavy(front)[0] ?? front[0], heavy(back)[0] ?? back[0]];
    const u = unit(sub(pos(mol, B), pos(mol, F)));
    const flat = v => sub(v, u.map(k => k * dot(v, u)));
    const e1 = unit(flat(sub(pos(mol, x[0]), pos(mol, F)))), e2 = cross(u, e1);
    const ang = (a, c) => { const v = sub(pos(mol, a), pos(mol, c)); return (Math.atan2(dot(v, e2), dot(v, e1)) * 180 / Math.PI + 360) % 360; };
    const phi0 = ang(x[1], B);
    const slot = (a, c, zero) => Math.round(((ang(a, c) - zero + 360) % 360) / 120) * 120 % 360; // 0, 120 or 240
    const side = new Set([B]);
    for (const a of side) for (const n of nbrs(mol, a)) if (n !== F) side.add(n); // the bond is not in a ring
    let levels = null, xx = null;
    const hasFlat = [...front, ...back].some(flatC);
    if ([...front, ...back].every(a => isH(a) || sp3C(a) || flatC(a))) {
      // xx: each non-H/non-H pair as [front atom, back atom, offset, gauche counts]; at slider angle phi their Newman
      // angle is phi + offset. Gauche counts only between two sp³ carbons.
      xx = heavy(front).flatMap(f => heavy(back).map(k => [f, k, slot(k, B, phi0) - slot(f, F, 0), sp3C(f) && sp3C(k)]));
      const key = [0, 60, 120, 180, 240, 300].map(phi => {
        const ecl = phi % 120 === 0;
        const n = xx.filter(([, , o, g]) => { const d = (phi + o + 720) % 360; return ecl ? d === 0 : g && (d === 60 || d === 300); }).length;
        return (ecl ? 100 : 0) + n;
      });
      const keys = [...new Set(key)].sort((p, q) => p - q);
      levels = key.map(k => keys.indexOf(k));
    }
    out.push({ atoms: [F, B], front, back, x, ref, phi0, side: [...side], levels, xx, flat: hasFlat });
  }
  return out;
}

// Steric (van der Waals) strain at slider angle phi: the non-H/non-H pairs of a ranked bond whose Newman angle is 65° or
// less, from gauche (60° ± 5°) down to eclipsed; a pair with a flat group only when eclipsed (5° or less), since its
// gauche form costs almost nothing. Returns [front atom, back atom, angle 0–180°] for each pair.
// Pairs with H are left out: in the textbook split (McMurry), CH3/CH3 gauche is steric strain, CH3/CH3 eclipsed is
// torsional + steric, and H/CH3 eclipsed is mostly torsional.
export function stericPairs(bond, phi) {
  return (bond.xx || []).map(([f, k, o, g]) => { const d = ((phi + o) % 360 + 360) % 360; return [f, k, Math.min(d, 360 - d), g]; })
    .filter(([, , a, g]) => a <= (g ? 65 : 5)).map(([f, k, a]) => [f, k, a]);
}

// The form at slider angle phi (the x-pair dihedral): staggered or eclipsed within 5°, otherwise in between (kind null).
// With reference groups, a staggered form is anti (180°) or gauche (60°, 300°).
export function conformer(phi, ref) {
  const k = Math.round(phi / 60);
  if (Math.abs(phi - 60 * k) > 5) return { kind: null, at: null, rel: null };
  const at = (k * 60) % 360, kind = k % 2 ? 'staggered' : 'eclipsed';
  return { kind, at, rel: ref && kind === 'staggered' ? (at === 180 ? 'anti' : 'gauche') : null };
}

// Each ring (OCL's ring set, up to 7 atoms): the angle at every ring atom between its two ring neighbors, taken from
// atomInfo's measured angles, with their mean and the mean ideal angle for the atoms' hybridization (null if unknown).
export function ringAngles(mol, atoms) {
  const rs = mol.getRingSet(), out = [];
  for (let r = 0; r < rs.getSize(); r++) {
    const ring = Array.from(rs.getRingAtoms(r)), n = ring.length; // plain array: OCL's cannot cross postMessage
    const angles = ring.map((a, i) => {
      const p = ring[(i + n - 1) % n], q = ring[(i + 1) % n];
      return atoms[a].angles.find(([s, t]) => (s === p && t === q) || (s === q && t === p))[2];
    });
    const hyb = ring.map(a => atoms[a].hyb ?? null);
    const mean = v => v.reduce((s, x) => s + x, 0) / v.length;
    out.push({ size: n, atoms: ring, aromatic: rs.isAromatic(r), hyb, angles, mean: mean(angles),
      ideal: hyb.every(Boolean) ? mean(hyb.map(h => IDEAL[h])) : null });
  }
  return out;
}

// ---- Shape tab: chairs (six-membered rings of sp³ atoms), axial/equatorial, ring flip ----
// Ring frame: center and unit normal (Newell: sum of cross products around the ring; atoms come in ring order).
function ringFrame(mol, ring) {
  const c = ring.map(a => pos(mol, a)).reduce((s, p) => s.map((v, i) => v + p[i] / ring.length), [0, 0, 0]);
  let n = [0, 0, 0];
  ring.forEach((a, i) => { n = n.map((v, k) => v + cross(sub(pos(mol, a), c), sub(pos(mol, ring[(i + 1) % ring.length]), c))[k]); });
  return { c, n: unit(n) };
}

// Signed dihedral a–b–c–d in degrees (−180 to 180).
function torsion(mol, a, b, c, d) {
  const [pa, pb, pc, pd] = [a, b, c, d].map(x => pos(mol, x)), b2 = sub(pc, pb);
  const n1 = cross(sub(pb, pa), b2), n2 = cross(b2, sub(pd, pc));
  return Math.atan2(dot(cross(n1, unit(b2)), n2), dot(n1, n2)) * 180 / Math.PI;
}

// A chair: the six ring torsions alternate in sign and each is at least 30° (ideal ±55°). A boat has two torsions
// near 0°, a twist-boat two neighbors of the same sign, a half-chair one near 0°.
export function isChair(mol, ring) {
  const t = ring.map((a, i) => torsion(mol, a, ring[(i + 1) % 6], ring[(i + 2) % 6], ring[(i + 3) % 6]));
  return t.every((x, i) => Math.abs(x) >= 30 && Math.sign(x) !== Math.sign(t[(i + 1) % 6]));
}

// Every bond from a ring atom to an atom outside the ring: axial when it runs within about 50° of the ring normal
// (ideal chair: axial 0°, equatorial about 70°), and which face it points to (up = along the normal).
export function chairSubs(mol, ring) {
  const { n } = ringFrame(mol, ring), inRing = new Set(ring);
  return ring.flatMap(a => nbrs(mol, a).filter(s => !inRing.has(s)).map(s => {
    const k = dot(unit(sub(pos(mol, s), pos(mol, a))), n);
    return { atom: a, sub: s, ax: Math.abs(k) > 0.65, up: k > 0, h: mol.getAtomicNo(s) === 1 };
  }));
}

// Why a chair is not flipped here (null = it can be): 'h' only H on the ring, so both chairs are the same;
// 'bridged'; 'trans' fused with the other ring on two equatorial bonds (locked, like trans-decalin); 'cis' fused
// with one axial bond (it can flip, but the other ring flips with it).
// shortcut: cis-fused rings (cis-decalin) are not flipped; fusion with a ring above 7 atoms (outside OCL's ring set)
// goes unseen, so the conformer search just finds no flipped chair there.
function flipBlock(mol, ring, subs) {
  if (subs.every(s => s.h)) return 'h';
  const rs = mol.getRingSet(), fused = [];
  for (let r = 0; r < rs.getSize(); r++) {
    const other = Array.from(rs.getRingAtoms(r)), shared = ring.filter(a => other.includes(a));
    if (shared.length === 6 || shared.length < 2) continue; // itself, or spiro: a spiro ring flips on its own
    const [p, q] = shared, i = ring.indexOf(p), j = ring.indexOf(q);
    if (shared.length > 2 || (Math.abs(i - j) !== 1 && Math.abs(i - j) !== 5)) return 'bridged';
    fused.push(subs.filter(s => (s.atom === p || s.atom === q) && other.includes(s.sub)).every(s => !s.ax));
  }
  return !fused.length ? null : fused.every(Boolean) ? 'trans' : 'cis';
}

// Which chair is more stable comes from MMFF94s+, so it is only told for groups whose axial/equatorial energy was
// checked against measured A-values in test.mjs. Alkyl (sp³ C with only C and H on it), OH and NH₂ come out on the
// measured side. F, Cl, Br and C≡N come out on the wrong side (axial lower); I comes out right, but as one family
// the halogens are left out. Returns null (a verdict can be given), {why: 'halogen' | 'nitrile' | 'unchecked'}.
function noVerdict(mol, subs) {
  for (const { atom, sub: s } of subs.filter(x => !x.h)) {
    const z = mol.getAtomicNo(s), nb = nbrs(mol, s), h = nb.filter(x => mol.getAtomicNo(x) === 1).length;
    if ([9, 17, 35, 53].includes(z)) return { why: 'halogen' };
    if (z === 6 && nb.length === 2 && piBonds(mol, s) === 2 && nb.some(x => mol.getAtomicNo(x) === 7)) return { why: 'nitrile' };
    const ok = (z === 6 && nb.length === 4 && nb.every(x => x === atom || [1, 6].includes(mol.getAtomicNo(x))))
      || (z === 8 && nb.length === 2 && h === 1) || (z === 7 && nb.length === 3 && h === 2);
    if (!ok) return { why: 'unchecked' };
  }
  return null;
}

// Chairs in the 3D model: six-membered rings (OCL's ring set) whose atoms are all sp³ and whose shape is a chair.
export function chairs(mol, atoms) {
  const rs = mol.getRingSet(), out = [];
  for (let r = 0; r < rs.getSize(); r++) {
    const ring = Array.from(rs.getRingAtoms(r));
    if (ring.length !== 6 || ring.some(a => atoms[a].hyb !== 'sp3') || !isChair(mol, ring)) continue;
    const subs = chairSubs(mol, ring);
    out.push({ atoms: ring, subs, block: flipBlock(mol, ring, subs), judge: noVerdict(mol, subs) });
  }
  return out;
}

// 1,3-diaxial pairs in a chair: each axial group other than H, with every axial atom two ring atoms away on the same
// face (the axial bonds of ring atoms 1, 3 and 5 all point to one face). Returns [group, partner] pairs, each once.
export function diaxial(chair) {
  const ring = chair.atoms, ax = chair.subs.filter(s => s.ax), out = [], seen = new Set();
  for (const s of ax.filter(x => !x.h)) {
    const i = ring.indexOf(s.atom);
    for (const t of ax.filter(x => x.up === s.up && [2, 4].includes((ring.indexOf(x.atom) - i + 6) % 6))) {
      const key = pairKey(s.sub, t.sub);
      if (!seen.has(key)) { seen.add(key); out.push([s.sub, t.sub]); }
    }
  }
  return out;
}

// The other chair of ring `ring`, from real conformers: OCL's conformer generator (a fresh start per seed) gives
// shapes that MMFF94s+ then relaxes fully; a run only counts when minimise() converged (returns 0). A conformer
// whose ring is a chair with every non-H group's axial/equatorial swapped is the flipped chair; with none swapped,
// the same chair. E = the lowest energy found for each chair (kJ/mol), the model included. `xyz` = the lowest flipped
// conformer whose groups keep their faces (a true flip, not its mirror image), turned to sit on the model's ring.
// Returns null when no flipped chair turns up.
// One conformer per seed (getOneConformerAsMolecule, about 30 ms; initializeConformers costs about 600 ms a seed).
// The ring shape comes from the generator's ring-fragment step, which changes with the seed. Tested: it gives both
// chairs for CH₃, OH, NH₂, halogens, C≡N and the dimethylcyclohexanes, but never the axial chair for ethyl or anything
// larger (tert-butyl, isopropyl, phenyl, OCH₃, COOH, vinyl, menthol, glucose), with any seed or strategy.
// shortcut: per chair the lowest of a few dozen minima, not a full search; a long flexible side chain can hide a
// lower minimum. Raise `seeds` if a verdict ever looks off.
export function flipChair(OCL, mol, ring, seeds = 40) {
  const heavy = s => s.filter(x => !x.h);
  const ref = heavy(chairSubs(mol, ring));
  const kind = m => { // 0 = same chair, 1 = flipped, -1 = neither
    if (!isChair(m, ring)) return -1;
    const now = heavy(chairSubs(m, ring));
    return now.every((s, i) => s.ax === ref[i].ax) ? 0 : now.every((s, i) => s.ax !== ref[i].ax) ? 1 : -1;
  };
  const E = [Infinity, Infinity];
  let best = null;
  const take = m => {
    m.ensureHelperArrays(OCL.Molecule.cHelperNeighbours);
    if (kind(m) < 0) return; // the ring shape rarely changes while minimising, so skip what cannot count
    const ff = new OCL.ForceFieldMMFF94(m, 'MMFF94s+');
    if (ff.minimise() !== 0) return;
    const k = kind(m), e = ff.getTotalEnergy() * 4.184;
    if (k < 0) return;
    E[k] = Math.min(E[k], e);
    if (k && heavy(chairSubs(m, ring)).every((s, i) => s.up === ref[i].up) && (!best || e < best.e)) best = { e, m };
  };
  take(mol.getCompactCopy());
  for (let seed = 1; seed <= seeds; seed++) {
    const m = new OCL.ConformerGenerator(seed).getOneConformerAsMolecule(mol.getCompactCopy());
    if (m) take(m);
  }
  if (!best) return null;
  return { E, dE: Number.isFinite(E[0]) ? E[1] - E[0] : null, xyz: alignRing(best.m, mol, ring), subs: chairSubs(best.m, ring) };
}

// Which chair is more stable from dE = (lowest flipped − lowest as-built) MMFF94s+ energy, kJ/mol: 'built', 'flipped',
// or 'same' when they differ by less than 1 kJ/mol.
export const moreStable = dE => (Math.abs(dE) < 1 ? 'same' : dE > 0 ? 'built' : 'flipped');

// Coordinates of every atom of `m`, moved and turned so its ring sits on the ring of `ref`: centers together,
// normals together, then the turn about the normal that best lines up the ring atoms. Only a rigid move.
function alignRing(m, ref, ring) {
  const f = ringFrame(m, ring), t = ringFrame(ref, ring);
  const rot = (v, axis, ang) => { // Rodrigues: turn v by ang about the unit axis
    const c = Math.cos(ang), s = Math.sin(ang), k = cross(axis, v), d = dot(axis, v);
    return v.map((x, i) => x * c + k[i] * s + axis[i] * d * (1 - c));
  };
  const ax = cross(f.n, t.n), sin = Math.hypot(...ax), a1 = Math.atan2(sin, dot(f.n, t.n));
  const axis = sin > 1e-9 ? unit(ax) : unit(cross(f.n, Math.abs(f.n[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0])); // normals (anti)parallel: any perpendicular axis
  const r1 = v => rot(v, axis, a1);
  let sc = 0, cc = 0;
  for (const a of ring) {
    const p = r1(sub(pos(m, a), f.c)), q = sub(pos(ref, a), t.c);
    sc += dot(cross(p, q), t.n); cc += dot(p, q);
  }
  const a2 = Math.atan2(sc, cc);
  return Array.from({ length: m.getAllAtoms() }, (_, i) => rot(r1(sub(pos(m, i), f.c)), t.n, a2).map((x, k) => x + t.c[k]));
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

// ---- Structure tab: polarity (bond dipoles and their vector sum) ----
// Pauling electronegativities.
export const EN = { H: 2.20, B: 2.04, C: 2.55, N: 3.04, O: 3.44, F: 3.98, Si: 1.90, P: 2.19, S: 2.58, Cl: 3.16, Se: 2.55, Br: 2.96, I: 2.66 };
// A bond counts as polar from a difference of 0.4: C–N (0.49) and C≡N must count (pyridine and acetonitrile are polar),
// C–H (0.35) must not.
export const POLAR_BOND = 0.4;
// The molecule is polar when the sum of the bond dipoles (each = its difference along the bond, δ+ → δ−) is at least
// this long. Symmetric molecules come out near 0 (below 0.02 in test.mjs), polar ones at 0.3 or more.
export const POLAR_NET = 0.1;
// A verdict only for molecules made of these: for C–I, C–S, C–P and S–H the difference understates the bond's polarity,
// and their lone pairs (left out here) matter more (CH₃I and dimethyl sulfide are polar).
const VERDICT_EL = new Set(['H', 'B', 'C', 'N', 'O', 'F', 'Cl', 'Br']);

// Returns { bonds: [[δ+ atom, δ− atom, difference]], net: [x, y, z], polar: true | false | null, why }.
// why (polar = null): 'ion' (net charge) or 'element' (an atom outside VERDICT_EL). Lone pairs are not counted.
export function polarity(mol) {
  const bonds = [];
  let net = [0, 0, 0], charge = 0, el = true;
  for (let a = 0; a < mol.getAllAtoms(); a++) {
    charge += mol.getAtomCharge(a);
    if (!VERDICT_EL.has(mol.getAtomLabel(a))) el = false;
  }
  for (let b = 0; b < mol.getAllBonds(); b++) {
    let [p, q] = [mol.getBondAtom(0, b), mol.getBondAtom(1, b)];
    const [ep, eq] = [EN[mol.getAtomLabel(p)], EN[mol.getAtomLabel(q)]];
    if (ep === undefined || eq === undefined) continue;
    const d = Math.round(Math.abs(ep - eq) * 100) / 100; // 3.44 − 3.04 is 0.3999… in floating point
    if (d < POLAR_BOND) continue;
    if (ep > eq) [p, q] = [q, p];
    bonds.push([p, q, d]);
    net = net.map((v, k) => v + d * unit(sub(pos(mol, q), pos(mol, p)))[k]);
  }
  const why = charge ? 'ion' : el ? null : 'element';
  return { bonds, net, polar: why ? null : Math.hypot(...net) >= POLAR_NET, why };
}

// ---- Structure tab: aromatic rings, π electrons and p orbitals ----
// Each ring atom gives its p orbital's electrons (Kekulé bonds from OCL): 1 for a double bond to an atom of the same
// aromatic ring system ('double'); 2 for a lone pair in the p orbital ('lp': no double bond, like pyrrole N–H, furan O,
// a ring C⁻); 0 for an empty p orbital ('empty': a ring C⁺). A ring with any other atom is left out.
// `planar` marks a 'double' atom whose lone pair lies in the ring plane (pyridine N), so it adds 0.
// Rings come from OCL's ring set (up to 7 atoms) where OCL calls them aromatic. Fused rings form a system, whose
// total counts each atom once. Returns { rings: [{ atoms, parts: [{ atom, e, why, planar }], pi, system }], systems: [{ rings, pi }] }.
export function aromaticRings(mol) {
  const rs = mol.getRingSet(), rings = [];
  for (let r = 0; r < rs.getSize(); r++) if (rs.isAromatic(r)) rings.push(Array.from(rs.getRingAtoms(r)));
  // ring systems: rings sharing an atom
  const sys = rings.map((_, i) => i), find = i => (sys[i] === i ? i : (sys[i] = find(sys[i])));
  rings.forEach((x, i) => rings.forEach((y, j) => { if (j > i && x.some(a => y.includes(a))) sys[find(j)] = find(i); }));
  const sysAtoms = i => new Set(rings.filter((_, j) => find(j) === find(i)).flat());
  const part = (a, inSys) => {
    const lp = lonePairs(mol, a);
    let dbl = false, exo = false;
    for (let k = 0; k < mol.getAllConnAtoms(a); k++) if (mol.getConnBondOrder(a, k) === 2) {
      if (inSys.has(mol.getConnAtom(a, k))) dbl = true; else exo = true;
    }
    if (dbl) return { atom: a, e: 1, why: 'double', planar: lp > 0 };
    if (exo) return null; // a C=O on the ring (2-pyridone): OCL does not call these rings aromatic anyway
    if (lp > 0) return { atom: a, e: 2, why: 'lp' };
    if (lp === 0 && mol.getAtomCharge(a) > 0) return { atom: a, e: 0, why: 'empty' };
    return null; // not a p-orbital atom we can count
  };
  const out = [], systems = new Map();
  rings.forEach((ring, i) => {
    const inSys = sysAtoms(i), parts = ring.map(a => part(a, inSys));
    if (parts.some(p => !p)) return;
    const s = find(i);
    if (!systems.has(s)) systems.set(s, { rings: [], atoms: inSys });
    systems.get(s).rings.push(out.length);
    out.push({ atoms: ring, parts, pi: parts.reduce((t, p) => t + p.e, 0), system: s });
  });
  const list = [...systems.entries()].map(([s, v]) => {
    const parts = [...v.atoms].map(a => part(a, v.atoms));
    return { key: s, rings: v.rings, pi: parts.every(Boolean) ? parts.reduce((t, p) => t + p.e, 0) : null };
  });
  out.forEach(r => { r.system = list.findIndex(x => x.key === r.system); });
  return { rings: out, systems: list.map(({ rings: r, pi }) => ({ rings: r, pi })) };
}

// 4n + 2 (Hückel): n, or null when the count is not 4n + 2.
export const huckel = pi => (pi >= 2 && (pi - 2) % 4 === 0 ? (pi - 2) / 4 : null);

// ---- Structure tab: chemically equivalent H (¹H NMR signals) ----
// The substitution test: replace one H with D and write OCL's canonical code (stereo included, read from the 3D model).
// Two H are equivalent when the two D compounds are the same (homotopic) or mirror images (enantiotopic); diastereomers
// (diastereotopic H) stay apart. The code ignores conformation, so fast motions (CH₃ rotation, ring flip) average out.
// Returns one class per H atom: { h: [atom], cls: [class per entry of h] }, classes numbered 0, 1, … by first appearance.
export function equivalentH(mol) {
  const h = [], keys = [];
  for (let a = 0; a < mol.getAllAtoms(); a++) {
    if (mol.getAtomicNo(a) !== 1) continue;
    const m = mol.getCompactCopy();
    m.setAtomMass(a, 2);
    const id = m.getIDCode();
    for (let i = 0; i < m.getAllAtoms(); i++) m.setAtomX(i, -m.getAtomX(i));
    h.push(a);
    keys.push([id, m.getIDCode()].sort().join(' '));
  }
  const seen = [];
  return { h, cls: keys.map(k => (seen.includes(k) ? seen : (seen.push(k), seen)).indexOf(k)) };
}
