# Section 6 auxiliary-slot geometry provenance

Status: **formal-structure / constructive paper-admissible witness**.  This file does not promote Stage 3 to paper-exact.

## Manuscript location

Section 6.2, especially the construction around Eqs. (6.9)-(6.14) and Lemma 6.1.  The manuscript assigns each slow label an auxiliary torus rectangle with a rational centre, uses a finite proper colouring of the interaction graph, bounds the covering-index difference of interacting bands, and then chooses one common positive rectangle radius `r0` so enlarged slots are quotient-injective and interacting labels have disjoint absolute lifted supports.

## Pinned official Lean cross-check

Repository: `openai/NavierStokesAndEuler`

Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant files and definitions:

- `NavierStokes/SlotColoring.lean`
  - `Palette = Fin 9 × ((Fin 3 → Fin 5) × Bool)` and `palette_card = 2250`;
  - `colorData`: dyadic level modulo 9, each integer grid coordinate modulo 5, and the sign bit;
  - `nativeGapBudget`, `nativeGap`, and theorem `nativeIndex_gap` for interacting levels separated by at most four.
- `NavierStokes/SlotGeometry.lean`
  - covering matrix `J=[[3,1],[1,5]]` and the bound `||J x||_∞ <= 6 ||x||_∞`;
  - `denominator(m,D)=(m+1)6^D`;
  - explicit rational centres `center(m,D,i)=((i+1)/denominator(m,D),0)`;
  - `centers_avoid_covering`, `exists_common_radius`, `orientedRectangle_subset`, `exists_oriented_slots`, and `colored_labels_have_slots`.
- `NavierStokes/PartitionedCovariance.lean`
  - `SlotSystem` and `exists_slotSystem`, which consume those centres/radius to obtain quotient injectivity and cross-label lifted-support disjointness.

The Lean development itself obtains the common radius through an open-neighbourhood existence argument.  It does not expose a canonical numerical value for that radius.

## Executable witness implemented here

`src/openai_ns_reconstruction/slot_geometry.py` makes the finite geometry reproducible without sampling a velocity field:

1. `palette_coordinates` implements the official modulo-9 / modulo-5 / sign colour data.
2. `palette_index` chooses a deterministic mixed-radix enumeration of the 2250 colours.
3. `native_gap(h)` implements the pinned `SlotColoring.nativeGap` formula and the pair certificate checks the actual landed `DyadicChart.covering_index` values.
4. `center_for_color` uses the exact rational centre formula from `SlotGeometry.center` with Python `Fraction` arithmetic.
5. The universal forbidden-centre torus gap is bounded below by `1/((2250+1)6^D)`: every centre and every integer-covering image has that common denominator, while the forbidden coincidences are nonzero modulo `Z^2` by the same positivity/injectivity mechanism formalized in `centers_avoid_covering`.
6. Instead of relying on the noncomputable neighbourhood radius, `ExplicitSlotSystem.standard_radius` chooses the smaller exact rational witness

   `R = min(1/(8*6^D), 1/(4*den*(6^D+1)))`.

   Hence `4 R 6^D < 1` gives a strict quotient-injectivity margin, while `2 R (6^D+1) < 1/den` preserves a strict cross-centre separation margin for every covering gap up to `D`.
7. For the manuscript axes `v_r=(1,-beta)`, `v_t=(beta,1)`, `beta=sqrt(2)-1`, both sup norms are one.  The exact `orientedRectangle_subset` denominator is therefore `1+1+1=3`, so the executable common coordinate radius is `r0=R/3`.
8. `pair_certificate` uses only exact rational/integer slot arithmetic plus the landed chart covering indices.  It does not inspect or fit base-field values.

Focused tests independently multiply the covering matrix, enumerate representative forbidden centre pairs in exact rational arithmetic, check the global radius inequalities, cross-check the chart covering-index gap, cover the two sign labels on one box and cross-level labels, and fail closed when the discrete interaction hypotheses are not certified.

## Deliberate non-equivalences and remaining gate

- The manuscript permits a finite proper colouring but does not prescribe a unique palette enumeration.  The pinned Lean uses `Fintype.equivFin Palette`; this repository uses an explicit mixed-radix enumeration.  The resulting centre assignment is therefore **paper-admissible**, not claimed definitionally identical to Lean's enumeration.
- The explicit rational `r0` is a constructive smaller witness satisfying the same separation estimates; it is not claimed equal to the radius hidden behind Lean's existential choice.
- `certify_interaction_colors` verifies the discrete consequences used by the colouring proof (`|ell-ell'|<=4`, and at equal level coordinate grid gaps `<=4`).  The remaining analytic/support obligation is to prove that the landed physical slow cutoffs/enlargements imply those hypotheses for every actually interacting pair, especially across different dyadic scalings.
- The paper-exact Proposition 5.5 background is still unavailable.  No arbitrary `TangentialBaseJetProvider`, phase sample, or slot certificate is promoted to a paper-exact base field.
- Uniform `LocalBaseBounds` / C1-C2 estimates and Eqs. (7.9)-(7.11) therefore remain blocked upstream.

Accordingly `paper_exact_velocity_available` and the Stage 3 status must remain false / `formal-structure`.
