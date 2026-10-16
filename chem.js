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
export function build3D(OCL, smiles) {
  // getFragments, not stripSmallFragments: the latter neutralizes charges (acetate would become acetic acid)
  const frags = OCL.Molecule.fromSmiles(smiles).getFragments();
  const mol = frags.reduce((a, b) => (b.getAtoms() > a.getAtoms() ? b : a));
  const stripped = frags.length > 1;
  const m3 = new OCL.ConformerGenerator(42).getOneConformerAsMolecule(mol);
  if (!m3) throw new Error('Could not build a 3D model for this structure.');
  let minimized = true;
  try { new OCL.ForceFieldMMFF94(m3, 'MMFF94s+').minimise(); }
  catch { minimized = false; } // atom types MMFF does not cover (e.g. metals): keep the raw conformer
  return { mol: m3, stripped, minimized };
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
