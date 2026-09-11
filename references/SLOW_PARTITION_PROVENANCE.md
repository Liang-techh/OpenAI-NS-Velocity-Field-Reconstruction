# Section 6 squared slow-partition provenance

Status: **formal-structure**. This increment constructs an explicit smooth squared partition satisfying the support and partition identities used in Section 6.2, but it does not upgrade Stage 3 to paper-exact and does not supply the missing Proposition 5.5 background or Lemma 6.1 rectangle coloring.

## Paper mapping

The implementation in `src/openai_ns_reconstruction/slow_partition.py` follows the construction immediately before Eqs. (6.8)-(6.9) of *Finite Time Blowup for Navier-Stokes* (Section 6.2, printed pages 64-65):

- choose a smooth squared partition `sum_ell chi_ell(q)^2 = 1`, with `chi_ell` supported where `q/Q in [1/2,2]`, `Q=2^-ell`;
- in each band choose a product squared partition `sum_a chi_{ell,a}(R,Z,T)^2 = 1` on mesh `S_*^-3`, with support extending by at most one mesh length on either side of its grid point in each coordinate;
- form the sign-independent slow cutoff `eta_gamma = chi_ell(q) chi_{ell,a}(R,Z,T)` in Eq. (6.9), so the `+` and `-` labels of one slow box duplicate the same cutoff.

The manuscript explicitly says both partitions can be constructed by normalizing translates of a smooth bump by the square root of their squared sum. It does **not** specify a unique seed bump or grid translate. The repository therefore chooses

`b(s)=exp(-1/(1-s^2))` for `|s|<1`, and `b(s)=0` otherwise,

and uses the zero-origin product grid. Those are paper-admissible implementation choices, not a claim of pointwise identity with an unpublished canonical cutoff.

For the dyadic partition the normalized lattice coordinate is `x=-log2(q)`, so `|x-ell|<=1` is exactly `q/Q in [1/2,2]`. For the slow product partition each coordinate is divided by `S_*^-3`. Because the seed bump is supported in `(-1,1)`, only two integer translates per one-dimensional coordinate can be nonzero; therefore the infinite squared sums reduce pointwise to two terms in the dyadic coordinate and at most eight terms in `(R,Z,T)`.

## Interface boundary

`ProductSlowCutoff.enclosure(label)` connects an actual product cutoff to the pre-existing `SlowBoxEnclosure` / `SlowLabel` interface and refuses mismatched band/grid labels. `slow_label_cutoff` also requires the caller's actual `DyadicChart`, preserving its `h` instead of silently recreating a default chart. The sign is intentionally absent from the cutoff value, matching Eq. (6.9).

This does not construct the centers `c_gamma` or the common radius `r0` in Lemma 6.1, does not prove the cross-label rectangle separation (6.13), and does not instantiate the base-field jet provider. Those remain independent blockers.

## Independent checks

`tests/test_slow_partition.py` verifies the squared partition identities by summing the independently enumerated local integer translates, checks the exact dyadic support endpoints, verifies the eight-term product localization, bridges a real product grid box to `SlowBoxEnclosure`, checks sign duplication, and exercises invalid/mismatched inputs fail-closed. These checks do not use a fitted wave or sampled background surrogate.

## Source pin

- Paper: OpenAI, *Finite Time Blowup for Navier-Stokes*, Section 6.2, Eqs. (6.8)-(6.9) and the squared-partition paragraph immediately preceding them.
- Official Lean repository pin retained by this project: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.
- No theorem in the pinned Lean snapshot is being claimed to select this particular bump or grid origin.
