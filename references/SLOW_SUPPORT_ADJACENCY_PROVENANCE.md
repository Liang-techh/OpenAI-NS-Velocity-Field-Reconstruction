# Section 6 slow-support to slot-adjacency provenance

Status: **formal-structure / geometric certificate**. This increment does not promote Stage 3 to paper-exact and does not instantiate the missing Proposition 5.5 background.

## Paper and official-Lean mapping

The slow cutoff before Eqs. (6.8)-(6.9) is a product of a dyadic factor supported on `q/Q in [1/2,2]` and a product factor supported within one `S_*^-3=ell^-6` mesh of an integer grid point in each normalized slow coordinate `(R,Z,T)`.

The project is pinned to `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`. In `NavierStokes/SlotColoring.lean`:

- `axisExponent D = [1/2,D,1]`;
- `spacing(a,n)=2^(-n a)/n^6` and `width D j n = spacing(axisExponent D j,n)`;
- `physicalBox` is the closed box of half-width `2*width` about the physical grid center;
- `Adj D L M` requires positive levels, distinct labels, level gap at most four, and a nonempty intersection of those physical boxes;
- `same_level_grid_gap` derives the integer coordinate gap `<=4` from that overlap.

For the manuscript similarity chart, `D=1/2-h` and `Q=2^-ell`. Therefore the physical images of one normalized product mesh are exactly

`Q^(1/2) ell^-6`, `Q^(1/2-h) ell^-6`, `Q ell^-6`,

which are the three official `SlotColoring.width` values. This is the cross-band scaling that must not be replaced by a same-scale comparison of normalized boxes.

## Executable bridge

`src/openai_ns_reconstruction/slow_support_adjacency.py` implements the following general certificate.

1. `physical_mesh_widths` specializes the official `spacing/width` formula to the landed chart exponents and rejects binary64 underflow/overflow instead of silently weakening the certificate.
2. `PhysicalSlowBox` converts one label's product-grid center to physical coordinates and measures support offsets by sending a *physical* point back through that label's own `DyadicChart`. Consequently two different dyadic levels are checked using their different physical scales.
3. The landed `ProductSlowCutoff` has support closure within one normalized mesh. The official Lean `physicalBox` has radius two meshes, so the entire landed support is contained globally and leaves one mesh of room for a separately justified fixed enlargement. The API permits an enlargement radius only in `(0,2]` and fails closed beyond the official envelope.
4. `certify_slow_support_adjacency` solves the shared physical similarity scale `q(z,tau)` once. Requiring the same physical point to lie in both dyadic support intervals `[Q/2,2Q]` yields the stronger exact integer consequence `|ell-ell'|<=2`, hence the `<=4` hypothesis required by `Adj`.
5. The same common physical point is converted through both charts and checked against each product-grid enclosure. It is therefore a witness of the nonempty `physicalBox` intersection. At equal level the implementation additionally checks the exact integer consequence `|a_j-b_j|<=4` before delegating to `slot_geometry.ExplicitSlotSystem.pair_certificate`.
6. The returned certificate therefore composes the landed slow-support geometry with the previously landed proper-coloring/rational-center/common-`r0` slot certificate without inspecting or fitting any base-field values.

## Independent regression checks

`tests/test_slow_support_adjacency.py` independently checks the official physical spacing formula, verifies one-mesh product support lies inside the two-mesh physical box, constructs a genuine cross-band common physical point and confirms each band uses its own scaling, exercises the sharp same-level grid-gap-four boundary for the two-mesh envelope, and verifies fail-closed behavior outside the dyadic support and radius hypotheses.

The cross-band test constructs the physical point first, maps that same point into both charts, and only then chooses the nearest product-grid indices. It does not compare normalized centers from different bands directly.

## Remaining gate

This closes the reusable **slow-support/enlargement -> `SlotColoring.Adj` -> explicit slot-separation** geometric adapter for any certified common support point inside the official two-mesh envelope. It does **not** supply:

- the paper-exact Proposition 5.5 background or an actual `TangentialBaseJetProvider`;
- theorem-selected `u_*` or any missing upstream profile data;
- uniform `LocalBaseBounds` / C1-C2 estimates on active boxes;
- the uniform Eqs. (7.9)-(7.11) phase/frame/damping theorem instance;
- stress-cone decomposition, amplitude ODEs, curl waves, mean corrections, or Section 9 residual iteration.

No caller-supplied field, sampled surrogate, generic wave, or finite point check is promoted to paper-exact data. `paper_exact_velocity_available` must remain `false`, and Stage 3 remains `formal-structure`.
