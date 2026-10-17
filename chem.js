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
  return mirrorID(mol) === mol.getIDCode();
}

// OCL's canonical code of the mirror image (x reflected; the code reads stereo from the 3D model).
function mirrorID(mol) {
  const m = mol.getCompactCopy();
  for (let a = 0; a < m.getAllAtoms(); a++) m.setAtomX(a, -m.getAtomX(a));
  return m.getIDCode();
}

// ---- Stereo tab: how two structures are related ----
// a, b = build3D results. In order: molecular formula; connectivity (canonical code with stereo removed: stripping it on
// a 3D copy does not work, since the code reads stereo back from the coordinates, so it goes through a copy without
// them); then stereo, from each 3D model's code and the code of b's mirror image. A meso compound is its own mirror
// image, so it comes out 'same' against itself. The code ignores conformation, so two shapes of one compound match.
// Returns { kind: 'not isomers' | 'constitutional' | 'open' | 'same' | 'enantiomers' | 'diastereomers', formulas, open }.
// 'open': a stereo unit is left open in a or b (open = [a, b]), checked only once the connectivity matches, since
// constitutional isomers are told apart without stereo (1-butanol and 2-butanol).
export function relation(OCL, a, b) {
  const formulas = [a, b].map(x => x.mol.getMolecularFormula().formula);
  if (formulas[0] !== formulas[1]) return { kind: 'not isomers', formulas };
  const flat = m => { const c = OCL.Molecule.fromIDCode(m.getIDCode()); c.stripStereoInformation(); return c.getIDCode(); };
  if (flat(a.mol) !== flat(b.mol)) return { kind: 'constitutional', formulas };
  const open = [a, b].map(stereoOpen);
  if (open.some(Boolean)) return { kind: 'open', formulas, open };
  const id = a.mol.getIDCode();
  return { kind: id === b.mol.getIDCode() ? 'same' : id === mirrorID(b.mol) ? 'enantiomers' : 'diastereomers', formulas };
}

// A stereocenter or stereo double bond that the name or SMILES left open (OCL parity 3), from build3D's `given`.
export const stereoOpen = ({ given }) => given.atoms.includes(3) || Object.values(given.bonds).includes(3);

// ---- Stereo tab: the mirror plane of a meso compound ----
// Stereocenters for meso: OCL's CIP centers, without pseudoasymmetric (r/s) ones, which sit on the plane.
export function cipCenters(OCL, mol) {
  mol.ensureHelperArrays(OCL.Molecule.cHelperCIP);
  const out = [];
  for (let a = 0; a < mol.getAtoms(); a++) if ([1, 2].includes(mol.getAtomCIPParity(a)) && !mol.isAtomParityPseudo(a)) out.push(a);
  return out;
}

// Meso: it matches its mirror image and has two or more stereocenters (cis-1,4-dimethylcyclohexane has none: achiral, not meso).
export const isMeso = (OCL, mol) => superimposable(mol) && cipCenters(OCL, mol).length >= 2;

// The substitution test on every atom: mark atom a with an isotope; a reflection that takes a to b makes "a marked" the
// mirror image of "b marked". partners[a] = every such b (all three H of a CH₃, since the code ignores conformation).
const MASS_NO = { 1: 1, 6: 12, 7: 14, 8: 16, 9: 19, 15: 31, 16: 32, 17: 35, 35: 79, 53: 127 };
function mirrorPartners(mol) {
  const n = mol.getAllAtoms(), key = [], mkey = [];
  for (let a = 0; a < n; a++) {
    const m = mol.getCompactCopy(), z = mol.getAtomicNo(a);
    m.setAtomMass(a, (MASS_NO[z] || 2 * z) + 1);
    key.push(m.getIDCode());
    for (let i = 0; i < n; i++) m.setAtomX(i, -m.getAtomX(i));
    mkey.push(m.getIDCode());
  }
  return mkey.map(k => key.flatMap((x, b) => (x === k ? [b] : [])));
}

const reflect = (p, c, n) => { const d = 2 * dot(sub(p, c), n); return p.map((x, k) => x - d * n[k]); };
const turnAbout = (P, from, to, atoms, deg) => { // positions with `atoms` turned by deg about the from → to bond
  const u = unit(sub(P[to], P[from])), th = deg * Math.PI / 180, c = Math.cos(th), s = Math.sin(th), o = P[to], Q = P.slice();
  for (const a of atoms) { const v = sub(P[a], o), k = cross(u, v), d = dot(u, v); Q[a] = v.map((x, i) => x * c + k[i] * s + u[i] * d * (1 - c) + o[i]); }
  return Q;
};

// How well the perpendicular bisector plane of atoms s and t reflects `atoms` onto their partners: the worst distance,
// the sum of squares, and whether each atom got its own partner.
function planeFit(P, part, s, t, atoms) {
  const c = P[s].map((x, k) => (x + P[t][k]) / 2), n = unit(sub(P[t], P[s])), used = new Set();
  let worst = 0, sum = 0;
  for (const a of atoms) {
    const r = reflect(P[a], c, n);
    let best = Infinity, b = -1;
    for (const j of part[a]) { const d = Math.hypot(...sub(P[j], r)); if (d < best) { best = d; b = j; } }
    worst = Math.max(worst, best); sum += best * best; used.add(b);
  }
  return { worst, sum, c, n, oneToOne: used.size === atoms.length };
}

// Turn single bonds outside rings until the plane between s and t fits, starting from positions P. The inner bonds of
// the chain from s to t (not its first or last bond, which only move side groups) go first, together on a 5° grid,
// scored on the chain atoms; then each other bond in turn over 0–355° in 5° steps and ±4° in 1° steps, heavy atoms
// first, then all atoms. Each bond turns its side away from s and t where it can, so the plane stays put.
// shortcut: three or more inner chain bonds are turned one at a time, which can miss; fine for the textbook meso set.
function turnToPlane(mol, P, part, s, t) {
  const side = (x, y) => { const S = new Set([x]); for (const a of S) for (const k of nbrs(mol, a)) if (!(a === x && k === y)) S.add(k); return S; };
  const bonds = [];
  for (let b = 0; b < mol.getAllBonds(); b++) {
    const [p, q] = [mol.getBondAtom(0, b), mol.getBondAtom(1, b)];
    if (mol.getBondOrder(b) !== 1 || mol.isRingBond(b) || mol.getAllConnAtoms(p) < 2 || mol.getAllConnAtoms(q) < 2) continue;
    const A = side(q, p), B = side(p, q), bad = S => [s, t].some(x => S.has(x) && x !== p && x !== q);
    const [from, to, S] = !bad(A) && (bad(B) || A.size <= B.size) ? [p, q, A] : !bad(B) ? [q, p, B] : A.size <= B.size ? [p, q, A] : [q, p, B];
    bonds.push({ from, to, atoms: [...S].filter(x => x !== to) });
  }
  const prev = { [s]: -1 }, queue = [s]; // shortest path from s to t
  for (const a of queue) for (const k of nbrs(mol, a)) if (!(k in prev)) { prev[k] = a; queue.push(k); }
  const path = []; for (let x = t; x !== -1; x = prev[x]) path.push(x);
  const onPath = new Set(path), inner = bonds.filter(({ from, to }) => onPath.has(from) && onPath.has(to) && ![s, t].includes(from) && ![s, t].includes(to));
  if (inner.length && inner.length <= 2) {
    let best = { e: Infinity, P };
    const grid = (Q, k) => {
      if (k === inner.length) { const e = planeFit(Q, part, s, t, path).sum; if (e < best.e) best = { e, P: Q }; return; }
      for (let d = 0; d < 360; d += 5) grid(turnAbout(Q, inner[k].from, inner[k].to, inner[k].atoms, d), k + 1);
    };
    grid(P, 0);
    P = best.P;
  }
  const all = [...P.keys()], heavy = all.filter(a => mol.getAtomicNo(a) !== 1);
  const rest = inner.length <= 2 ? bonds.filter(b => !inner.includes(b)) : bonds;
  for (const atoms of [heavy, heavy, all, all, all, all]) for (const { from, to, atoms: S } of rest) {
    const score = d => planeFit(turnAbout(P, from, to, S, d), part, s, t, atoms).sum;
    let best = { e: score(0), d: 0 };
    for (let d = 5; d < 360; d += 5) { const e = score(d); if (e < best.e - 1e-9) best = { e, d }; }
    for (let d = best.d - 4; d <= best.d + 4; d++) { const e = score(d); if (e < best.e - 1e-9) best = { e, d }; }
    if (best.d) P = turnAbout(P, from, to, S, best.d);
  }
  return P;
}

// Rigid move of positions P onto Q (least squares over all atoms; Horn's quaternion method, power iteration).
function fitOnto(P, Q) {
  const mean = X => [0, 1, 2].map(k => X.reduce((t, p) => t + p[k], 0) / X.length), cp = mean(P), cq = mean(Q);
  const S = [0, 1, 2].map(i => [0, 1, 2].map(j => P.reduce((t, p, a) => t + (p[i] - cp[i]) * (Q[a][j] - cq[j]), 0)));
  const [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]] = S;
  const N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx], [yz - zy, xx - yy - zz, xy + yx, zx + xz],
    [zx - xz, xy + yx, -xx + yy - zz, yz + zy], [xy - yx, zx + xz, yz + zy, -xx - yy + zz]];
  const shift = N.flat().reduce((t, v) => t + Math.abs(v), 0) + 1; // makes the largest eigenvalue the dominant one
  let q = [1, 0.1, 0.1, 0.1];
  for (let it = 0; it < 500; it++) q = unit(N.map((row, i) => row.reduce((t, v, j) => t + v * q[j], 0) + shift * q[i]));
  const [w, ...u] = q;
  return P.map(p => { const v = sub(p, cp), c1 = cross(u, v), c2 = cross(u, c1); return v.map((x, k) => x + 2 * w * c1[k] + 2 * c2[k] + cq[k]); });
}

// The mirror plane of a meso compound, checked on real coordinates: each atom, reflected through the plane, lands within
// 0.3 Å of a partner atom, one partner each, and no two atoms three or more bonds apart come closer than 1.5 Å.
// The plane swaps stereocenters (R ↔ S), so it is the perpendicular bisector of such a pair. Tried in order: the model
// as built ('model'); the model turned about single bonds ('turned', as on the Shape tab); conformers from OCL's
// generator (seeds 1–`seeds`, MMFF94s+, converged runs only, same stereo), turned the same way and moved onto the model
// ('conformer'). Rings never turn, so a ring keeps its shape: no chair of cis-1,2-dimethylcyclohexane has the plane.
// Returns { xyz, c, n, how, worst } or null.
export function mirrorPlane(OCL, mol, seeds = 20) {
  if (!isMeso(OCL, mol)) return null;
  const part = mirrorPartners(mol), centers = cipCenters(OCL, mol), all = Array.from({ length: mol.getAllAtoms() }, (_, a) => a);
  const pairs = centers.flatMap(s => part[s].filter(t => t > s && centers.includes(t)).map(t => [s, t]));
  const apart = all.map(a => { const near = new Set([a, ...nbrs(mol, a)]); for (const x of nbrs(mol, a)) for (const y of nbrs(mol, x)) near.add(y); return near; });
  const roomy = P => all.every(a => all.every(b => b <= a || apart[a].has(b) || Math.hypot(...sub(P[a], P[b])) >= 1.5));
  const check = P => {
    for (const [s, t] of pairs) {
      const f = planeFit(P, part, s, t, all);
      if (f.worst <= 0.3 && f.oneToOne && roomy(P)) return { xyz: P, c: f.c, n: f.n, worst: f.worst };
    }
    return null;
  };
  const try3 = (P, how) => {
    const now = check(P);
    if (now) return { ...now, how };
    for (const [s, t] of pairs) { const r = check(turnToPlane(mol, P, part, s, t)); if (r) return { ...r, how: how === 'model' ? 'turned' : how }; }
    return null;
  };
  const home = all.map(a => pos(mol, a)), id = mol.getIDCode();
  const found = try3(home, 'model');
  if (found) return found;
  for (let seed = 1; seed <= seeds; seed++) {
    const m = new OCL.ConformerGenerator(seed).getOneConformerAsMolecule(mol.getCompactCopy());
    if (!m) continue;
    m.ensureHelperArrays(OCL.Molecule.cHelperNeighbours);
    try { if (new OCL.ForceFieldMMFF94(m, 'MMFF94s+').minimise() !== 0) continue; } catch { continue; }
    if (m.getIDCode() !== id) continue;
    const r = try3(fitOnto(all.map(a => pos(m, a)), home), 'conformer');
    if (r) { const f = fitOnto(r.xyz, home), g = check(f); if (g) return { ...g, how: 'conformer' }; } // move the found shape onto the model, then check again
  }
  return null;
}

// ---- Stereo tab: Fischer projection ----
// The longest chain of carbons outside rings that holds every chain stereocenter, drawn top to bottom. Each stereocenter
// is a cross: its chain neighbors on the vertical line (pointing back, away from the viewer), its other two groups on the
// horizontal line (pointing out, toward the viewer). Which group goes right comes from the 3D model: with up U, down D and
// group A around center c, det(U − c, D − c, A − c) < 0 puts A on the right (in the ideal Fischer frame, U = (0, 1, −1),
// D = (0, −1, −1), right = (1, 0, 1), and the determinant is −2). That is the configuration, so any conformation gives it.
// Top end: the more oxidized one (more bonds to O); a tie goes to the end IUPAC numbers 1 (lowest locants for the
// principal group, then multiple bonds, then all substituents).
// shortcut: no IUPAC alphabetical tie-break, and a ring group is written as its formula (C₆H₅); add if a class case needs it.
// Returns { why: 'none' | 'ring' | 'chain' | 'end' } or { rows, dl, mirror, ringLeft, open }:
// rows top to bottom: { atom, label } or { atom, cross: true, up, down, left, right: { atom, label }, rs, row: index in `rows` of stereo() };
// dl = { dl: 'D' | 'L', kind: 'sugar' | 'amino', at: row index } for open-chain sugars (the lowest stereocenter) and α-amino
// acids (the α-carbon: L-threonine is L by C2, not C3); mirror = the row (or half row) the meso mirror plane crosses.
export function fischer(OCL, mol, rows) {
  mol.ensureHelperArrays(OCL.Molecule.cHelperCIP);
  const z = a => mol.getAtomicNo(a), heavy = a => nbrs(mol, a).filter(x => z(x) !== 1), Hs = a => nbrs(mol, a).filter(x => z(x) === 1).length;
  const bond = (a, b) => mol.getBondOrder(mol.getBond(a, b));
  const stereoAtom = a => a < mol.getAtoms() && [1, 2].includes(mol.getAtomParity(a)) && nbrs(mol, a).length === 4;
  const centers = Array.from({ length: mol.getAtoms() }, (_, a) => a).filter(stereoAtom);
  if (!centers.length) return { why: 'none' };
  const chainC = a => z(a) === 6 && !mol.isRingAtom(a);
  const xs = centers.filter(chainC); // the crosses
  if (!xs.length) return { why: centers.every(a => mol.isRingAtom(a)) ? 'ring' : 'chain' };
  // longest path of chain carbons holding every cross (they form a forest, so a path is unique per pair of ends)
  const pathFrom = u => { const prev = { [u]: -1 }, q = [u]; for (const a of q) for (const k of heavy(a)) if (chainC(k) && !(k in prev)) { prev[k] = a; q.push(k); } return prev; };
  const subs = p => p.reduce((t, a) => t + heavy(a).filter(x => !p.includes(x)).filter(x => z(x) !== 6).length, 0);
  let best = null;
  for (const u of Object.keys(pathFrom(xs[0])).map(Number)) {
    const prev = pathFrom(u);
    for (const v of Object.keys(prev).map(Number)) {
      const p = []; for (let x = v; x !== -1; x = prev[x]) p.push(x);
      if (!xs.every(c => p.includes(c))) continue;
      if (!best || p.length > best.length || (p.length === best.length && subs(p) > subs(best))) best = p;
    }
  }
  if (!best) return { why: 'chain' };
  let d = best;
  for (const e of [0, 1]) { // a cross at an end of the chain takes its ring-carbon neighbor as the next vertical atom
    const end = e ? d[d.length - 1] : d[0];
    if (!xs.includes(end)) continue;
    const next = heavy(end).find(x => z(x) === 6 && !d.includes(x));
    if (next === undefined) return { why: 'end' };
    d = e ? [...d, next] : [next, ...d];
  }
  // top end
  const oBonds = a => nbrs(mol, a).filter(x => z(x) === 8).reduce((t, x) => t + bond(a, x), 0);
  const rank = a => { // principal characteristic group on chain carbon a (lower ranks first): acid, nitrile, aldehyde, ketone, alcohol, thiol, amine
    const nb = heavy(a), dblO = nb.some(x => z(x) === 8 && bond(a, x) === 2), sglO = nb.some(x => z(x) === 8 && bond(a, x) === 1);
    if (dblO && sglO) return 1;
    if (nb.some(x => z(x) === 7 && bond(a, x) === 3)) return 3;
    if (dblO) return Hs(a) ? 4 : 5;
    if (nb.some(x => z(x) === 8 && Hs(x))) return 6;
    if (nb.some(x => z(x) === 16 && Hs(x))) return 7;
    if (nb.some(x => z(x) === 7)) return 8;
    return 99;
  };
  const locants = p => {
    const top = Math.min(...p.map(rank));
    return [p.flatMap((a, i) => (top < 99 && rank(a) === top ? [i + 1] : [])),
      p.flatMap((a, i) => (i + 1 < p.length && bond(a, p[i + 1]) > 1 ? [i + 1] : [])),
      p.flatMap((a, i) => heavy(a).filter(x => !p.includes(x)).map(() => i + 1))];
  };
  const lower = (x, y) => { for (let k = 0; k < x.length; k++) for (let i = 0; i < Math.max(x[k].length, y[k].length); i++) if (x[k][i] !== y[k][i]) return (x[k][i] ?? Infinity) < (y[k][i] ?? Infinity) ? -1 : 1; return 0; };
  const r = [...d].reverse();
  if (oBonds(d.at(-1)) > oBonds(d[0]) || (oBonds(d.at(-1)) === oBonds(d[0]) && lower(locants(r), locants(d)) < 0)) d = r;

  // rows: crosses, plus chain atoms between them; a pure alkyl tail (only C and H, single bonds, no ring) is one label
  const label = a => condensed(mol, a, d);
  const first = d.findIndex(a => xs.includes(a)), last = d.length - 1 - [...d].reverse().findIndex(a => xs.includes(a));
  const alkyl = (a, from) => { const S = new Set([a]); for (const x of S) for (const k of nbrs(mol, x)) if (k !== from && !S.has(k)) S.add(k);
    return [...S].every(x => [1, 6].includes(z(x)) && !mol.isRingAtom(x) && nbrs(mol, x).every(k => bond(x, k) === 1)); };
  const out = [], det = (u, v, w) => dot(u, cross(v, w));
  const push = i => {
    const a = d[i];
    if (!xs.includes(a)) { out.push({ atom: a, label: i === 0 || i === d.length - 1 ? condensed(mol, a, [d[i === 0 ? 1 : i - 1]]) : rowLabel(mol, a, d) }); return; }
    const [U, D] = [d[i - 1], d[i + 1]], [A, B] = nbrs(mol, a).filter(x => x !== U && x !== D), c = pos(mol, a);
    const right = det(sub(pos(mol, U), c), sub(pos(mol, D), c), sub(pos(mol, A), c)) < 0 ? A : B, left = right === A ? B : A;
    const row = rows.findIndex(x => x.kind === 'center' && x.atom === a);
    out.push({ atom: a, cross: true, up: U, down: D, left: { atom: left, label: label(left) }, right: { atom: right, label: label(right) }, rs: rows[row]?.rs ?? null, row });
  };
  if (first > 1 && alkyl(d[first - 1], d[first])) out.push({ atom: d[first - 1], label: condensed(mol, d[first - 1], [d[first]]) });
  else for (let i = 0; i < first; i++) push(i);
  for (let i = first; i <= last; i++) push(i);
  if (last < d.length - 2 && alkyl(d[last + 1], d[last])) out.push({ atom: d[last + 1], label: condensed(mol, d[last + 1], [d[last]]) });
  else for (let i = last + 1; i < d.length; i++) push(i);

  // D/L
  const OH = a => heavy(a).filter(x => z(x) === 8 && bond(a, x) === 1 && Hs(x) === 1);
  const isCHO = a => heavy(a).length === 2 && Hs(a) === 1 && heavy(a).some(x => z(x) === 8 && bond(a, x) === 2);
  const isCH2OH = a => heavy(a).length === 2 && Hs(a) === 2 && OH(a).length === 1;
  const isCO = a => heavy(a).length === 3 && heavy(a).some(x => z(x) === 8 && bond(a, x) === 2);
  const top = isCHO(d[0]) ? 1 : isCH2OH(d[0]) && isCO(d[1]) ? 2 : 0;
  const crosses = out.flatMap((x, i) => (x.cross ? [i] : []));
  let dl = null;
  if (d.length >= 3 && top && isCH2OH(d.at(-1)) && d.slice(top, -1).every(a => OH(a).length === 1 && Hs(a) === 1 && heavy(a).length === 3) && crosses.length) {
    const k = crosses.at(-1), x = out[k];
    dl = { kind: 'sugar', at: k, dl: OH(x.atom).includes(x.right.atom) ? 'D' : 'L' };
  }
  const carboxyl = a => heavy(a).length === 3 && heavy(a).some(x => z(x) === 8 && bond(a, x) === 2) && heavy(a).some(x => z(x) === 8 && bond(a, x) === 1 && heavy(x).length === 1);
  if (!dl && carboxyl(d[0]) && out[1]?.cross && out[1].atom === d[1]) {
    const N = heavy(d[1]).find(x => z(x) === 7 && heavy(x).length === 1);
    if (N !== undefined) dl = { kind: 'amino', at: 1, dl: out[1].right.atom === N ? 'D' : 'L' };
  }
  // meso look: the rows read the same from the bottom up (a meso compound drawn this way has its plane across the middle)
  const same = (x, y) => (x.cross ? y.cross && x.left.label === y.left.label && x.right.label === y.right.label : !y.cross && x.label === y.label);
  const mirror = isMeso(OCL, mol) && out.every((x, i) => same(x, out[out.length - 1 - i])) ? (out.length - 1) / 2 : null;
  return { rows: out, dl, mirror, ringLeft: centers.some(a => mol.isRingAtom(a)), open: out.some(x => x.cross && x.row >= 0 && !rows[x.row].given) };
}
// A group as a condensed formula, from atom a away from the atoms in `away`: CH₃, CH₂OH, CHO, COOH, CH(CH₃)₂, OCH₃, NH₂.
// A carbon writes its H, then a double-bonded O, then its other groups (all but the largest in parentheses, alike ones
// counted: (CH₃)₂). A group with a ring atom is written as its formula (C₆H₅). Digits stay digits; the page sets them low.
function condensed(mol, a, away) {
  const z = x => mol.getAtomicNo(x), out = new Set(away);
  const S = new Set([a]); for (const x of S) for (const k of nbrs(mol, x)) if (!out.has(k) && !S.has(k)) S.add(k);
  if ([...S].some(x => mol.isRingAtom(x))) {
    const n = {}; for (const x of S) n[mol.getAtomLabel(x)] = (n[mol.getAtomLabel(x)] || 0) + 1;
    const els = Object.keys(n).sort((p, q) => (p === 'C' ? -1 : q === 'C' ? 1 : p === 'H' ? -1 : q === 'H' ? 1 : p.localeCompare(q)));
    return els.map(e => e + (n[e] > 1 ? n[e] : '')).join('');
  }
  const write = (x, from) => {
    const nb = nbrs(mol, x).filter(k => k !== from && !out.has(k)), h = nb.filter(k => z(k) === 1).length;
    const q = mol.getAtomCharge(x), bo = k => mol.getBondOrder(mol.getBond(x, k));
    const oxo = nb.filter(k => z(k) === 8 && bo(k) === 2 && nbrs(mol, k).length === 1);
    const kids = nb.filter(k => z(k) !== 1 && !oxo.includes(k)).map(k => [(bo(k) === 2 ? '=' : bo(k) === 3 ? '≡' : '') + write(k, x), size(k, x)]).sort((p, r) => p[1] - r[1]);
    let s = mol.getAtomLabel(x) + (h ? 'H' + (h > 1 ? h : '') : '') + (q ? (Math.abs(q) > 1 ? Math.abs(q) : '') + (q > 0 ? '⁺' : '⁻') : '') + 'O'.repeat(oxo.length);
    for (let i = 0; i < kids.length; i++) {
      const t = kids[i][0], n = kids.filter(k => k[0] === t).length;
      if (n > 1) { s += `(${t})${n}`; i += n - 1; } else s += i < kids.length - 1 ? `(${t})` : t;
    }
    return s;
  };
  const size = (x, from) => { const T = new Set([x]); for (const y of T) for (const k of nbrs(mol, y)) if (k !== from && !out.has(k) && !T.has(k) && z(k) !== 1) T.add(k); return T.size; };
  return write(a, -1);
}

// A chain atom between crosses: C, its H, =O, then its side groups in parentheses: CH₂, C=O, CH(CH₃).
function rowLabel(mol, a, chain) {
  const nb = nbrs(mol, a).filter(k => !chain.includes(k)), z = x => mol.getAtomicNo(x);
  const oxo = nb.some(k => z(k) === 8 && mol.getBondOrder(mol.getBond(a, k)) === 2);
  const h = nb.filter(k => z(k) === 1).length, side = nb.filter(k => z(k) !== 1 && !(z(k) === 8 && mol.getBondOrder(mol.getBond(a, k)) === 2)).map(k => condensed(mol, k, [a]));
  const kinds = [...new Set(side)];
  return 'C' + (h ? 'H' + (h > 1 ? h : '') : '') + (oxo ? '=O' : '') + kinds.map(t => { const n = side.filter(x => x === t).length; return `(${t})${n > 1 ? n : ''}`; }).join('');
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
