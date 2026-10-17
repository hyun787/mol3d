# 3c progress (Reactions tab)

Working notes for the 3c branch. Deleted when 3c is merged.

## Done
- Group 1, backside attack (SN2): C–X rows (sp³ C with Cl, Br, I; equivalent ones once, mirror images included), textbook rate rule (`sn2Rate`), backside cone (60°, 3.0 Å from C), OH⁻ arrow, inverted alcohol as model 1 (placed on the inverted reactant), R/S before and after, guess first for the product's R/S. Reasons when there is no C–X (C–F, C–OH, vinyl/aryl halide).
- Group 2, flat intermediate (SN1, E1): textbook rule (`sn1Rate`: 3° good, 2° slow, 1°/methyl none, allylic/benzylic by resonance), carbocation built from the model with X off and C⁺ (converged MMFF94s+ only), flat check (angle sum ≥ 355°), lobes and two water arrows, water on each face → two alcohols as models 1 and 2, relation() → one product / enantiomers (racemic) / diastereomers, hydride/methyl shift note (graph rule: a neighbor would become a more substituted cation), allylic far-end note, radical note, guess first for the relation.
- Group 3, anti-periplanar (E2): β-H rows (H on one β carbon giving the same alkene share a row), open chain: the β side turned rigidly to H–C–C–X 180° (morph, Newman-like view tilted 25°), ring: trans-diaxial by configuration (faces), "ring flip puts X axial" with the Shape tab's flip when available, alkene per H from the anti shape (E/Z read by OCL from 3D), list of alkenes that can form, Zaitsev tag (only when alkenes differ in substitution), note when the Zaitsev alkene cannot form, Hofmann note, guess first (which alkene, or none).
- Shape tab fix found on the way: flipChair's conformer swaps the two H of some ring CH₂ (same shape, H identity lost); each H is now put back on its own face. test.mjs checks it.

## Next
- Group 4 (faces of C=C), guess-first pass, then README, merge, deploy, live check.

## Decisions taken (defaults, no user input)
- E2 rows: a ring bond in a ring that is not a chair in the model gets no verdict (the brief's rule is for chairs).
- Rate verdicts only when the carbon's other groups are C or H (1-bromo-1-chloroethane: no verdict, product still shown).
- C–X rows dedupe homotopic and enantiotopic C–X (substitution test with the mirror image), as for equivalent H.
- SN2 umbrella animation (optional in the brief) skipped for now: reactant and product have different atoms (Br vs O, the new O–H), so a straight morph needs an atom map and a made-up path.
- MMFF94s+ builds ring 3° carbocations at 354.5–354.7° (1-methylcyclohexyl, also from the 1,3-dimethyl case), below the brief's 355°: those cations are not drawn (text says so, with the number), but their two faces still give the products, since an O placed above or below sets the configuration either way. Open-chain 3° (356–359°), benzyl and allyl (360°) pass.

## Library behavior checked
- MMFF94s+ handles carbocations: tert-butyl, 3-methylhexan-3-yl, benzyl, allyl, isopropyl cations all converge (minimise() 0), angle sums 356–360°.
- ConformerGenerator works on an edited 3D molecule with explicit H: heavy-atom order kept, new H appended, stereo read from the edited coordinates (getIDCode of the edit = of the build).
- 3Dmol removeModel pops trailing ids, so removing the highest model first lets addModel reuse 1, 2.
- OCL toIsomericSmiles can write conflicting / and \ when a C=C carbon has a group on each side (Br/C(\c1ccccc1)=C/c1ccccc1); RDKit reads it without E/Z. Products are posted as molfiles, and test.mjs has RDKit read the molfile (stereo from 3D, removeHs) instead of OCL's SMILES.
- RDKit.js: get_mol(molfile, '{"removeHs":true}') gives the canonical SMILES with stereo from the 3D coordinates.
- 3Dmol animated zoomTo adds per-step deltas to the view: two at once, or one during our own animateView, overshoot. Instant zoomTo finds the target view; animateView moves there.

## Tabs width
- Five tabs at 280 px with 4 px side padding fit with the fallback font in the harness (scrollWidth 278). Recheck with IBM Plex Sans on the live site; rename to "React" if it overflows.
