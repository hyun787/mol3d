# 3c progress (Reactions tab)

Working notes for the 3c branch. Deleted when 3c is merged.

## Done
- Group 1, backside attack (SN2): C–X rows (sp³ C with Cl, Br, I; equivalent ones once, mirror images included), textbook rate rule (`sn2Rate`), backside cone (60°, 3.0 Å from C), OH⁻ arrow, inverted alcohol as model 1 (placed on the inverted reactant), R/S before and after, guess first for the product's R/S. Reasons when there is no C–X (C–F, C–OH, vinyl/aryl halide).

## Next
- Group 2 (SN1 flat carbocation), group 3 (E2), group 4 (faces of C=C), then README, merge, deploy, live check.

## Decisions taken (defaults, no user input)
- Rate verdicts only when the carbon's other groups are C or H (1-bromo-1-chloroethane: no verdict, product still shown).
- C–X rows dedupe homotopic and enantiotopic C–X (substitution test with the mirror image), as for equivalent H.
- SN2 umbrella animation (optional in the brief) skipped for now: reactant and product have different atoms (Br vs O, the new O–H), so a straight morph needs an atom map and a made-up path.

## Library behavior checked
- MMFF94s+ handles carbocations: tert-butyl, 3-methylhexan-3-yl, benzyl, allyl, isopropyl cations all converge (minimise() 0), angle sums 356–360°.
- ConformerGenerator works on an edited 3D molecule with explicit H: heavy-atom order kept, new H appended, stereo read from the edited coordinates (getIDCode of the edit = of the build).
- 3Dmol removeModel pops trailing ids, so removing the highest model first lets addModel reuse 1, 2.
- 3Dmol animated zoomTo adds per-step deltas to the view: two at once, or one during our own animateView, overshoot. Instant zoomTo finds the target view; animateView moves there.

## Tabs width
- Five tabs at 280 px with 4 px side padding fit with the fallback font in the harness (scrollWidth 278). Recheck with IBM Plex Sans on the live site; rename to "React" if it overflows.
