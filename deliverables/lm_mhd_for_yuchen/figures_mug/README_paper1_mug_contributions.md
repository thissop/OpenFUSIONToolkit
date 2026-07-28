# MUG contributions to paper #1 (finite-c QS-validity map)

Figures + text points MUG provides; box 1 integrates into paper1.tex.

## Figures (rendered, in this dir)
1. **`mug_loads_secondcoding.png`** — MUG-full vs COMSOL-full disruption loads at the headline
   (tau_q=1ms, c=1, Ha=2646, packing=800 run): intEJ total 0.1%, Iw 1.2%, EJw 0.9%, EJf 0.4%.
   Supports the sec:map "independent full-induction second code" paragraph (already in the tex).
2. **`mug_thinwall_boundary.png`** — MUG-thin vs COMSOL-resolved-wall Iw discrepancy vs tau_q:
   0.7% (10ms) -> 5.2% (3ms) -> 13.6% (1ms), monotonic across the delta_w/tw~1 crossover at 1.76ms.
   MUG's unique result: thin-wall BCs become unsafe precisely in the fast-disruption regime.
   Candidate for a paper-#1 subsection OR small-paper #4.

## Text points MUG certifies (for the map's claims)
- Loads two-code-robust to 0.1-1.2% (blind pre-registered), mesh-robust across 256^2/384^2.
- Full-vs-QS crossover (sign change of full-qs Umax at c~10) independently seen from the full side
  -> corroborates box 1's error-cancellation mechanism for the non-monotonic-c dip.
- MUG cannot two-code the dip itself (no QS mode; Umax odd-even artifact; loads monotonic) -> dip is
  well-resolved single-code + mechanism (c_dip~Ha^+1/2), NOT two-code, and paper should say so.
