# Section 7 exact large-band phase-scale gate provenance

Status: **formal-structure**. This artifact certifies only the scalar asymptotic band hypotheses feeding the pinned phase-error theorem. It does not construct the paper-exact base field, a complete slow-box estimate, or an oscillatory velocity wave.

## Pinned sources

Official Lean repository: `openai/NavierStokesAndEuler`

Pinned commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Relevant modules and statements:

- `NavierStokes/ChartScales.lean`
  - `Q n = SlotColoring.dyadicQ n`, the dyadic manuscript scale.
  - `S n = n^2`.
  - `epsilon h n = Q n ^ h`.
  - `carrier h n = Scaling.carrierFrequency (epsilon h n)`.
- `NavierStokes/Scaling.lean`
  - `carrierFrequency epsilon = ceil(epsilon^(-1/2))`.
  - `carrier_frequency_sqrt_bounds`, in particular the ceiling-rounded carrier is at least the continuum inverse-square-root frequency.
- `NavierStokes/PhaseEstimates.lean`
  - `phaseError S epsilon k = 1/S + S epsilon^2 + S/k + epsilon S`.
  - `phaseError_le_four_div`: if `S^2 epsilon^2 <= 1`, `S^2/k <= 1`, and `epsilon S^2 <= 1`, then `phaseError <= 4/S`.
  - `eventually_band_conditions` and `eventually_phaseError_le`: those hypotheses eventually hold on the actual dyadic scale and actual ceiling-rounded carrier.

Paper locations: Section 6 dyadic chart scale and Section 7.1 phase construction / Eqs. (7.9)–(7.11) estimate layer.

## Implemented mapping

Implemented in `src/openai_ns_reconstruction/phase_large_band_scale.py`.

For a positive rational `h=a/b` and integer band index `ell`, the module keeps the manuscript choices symbolically:

- `S_* = ell^2`;
- `epsilon = 2^(-ell h)` through the exact stored quantity `log2(epsilon) = -ell h`;
- `k = ceil(epsilon^(-1/2))` through the rigorous lower bound `k >= epsilon^(-1/2)`.

It never forms `Q=2^-ell` as binary64. The three scalar gates are reduced to integer comparisons:

1. `S_*^2 epsilon^2 <= 1` is equivalent to
   `ell^(4b) <= 2^(2 ell a)`.
2. Since `k >= epsilon^(-1/2)`, the sufficient condition
   `S_*^2 sqrt(epsilon) <= 1` implies `S_*^2/k <= 1`; after clearing the rational exponent this is
   `ell^(8b) <= 2^(ell a)`.
3. `epsilon S_*^2 <= 1` is equivalent to
   `ell^(4b) <= 2^(ell a)`.

The carrier comparison in item 2 is intentionally recorded as **sufficient, not necessary**: the integer ceiling can make `S_*^2/k` smaller than the continuum upper bound. No claim is made that the first passing band of this sufficient gate is the first band satisfying the exact ceiling-rounded theorem hypothesis.

For the repository default `h=1/200`, the independent regression checks that this sufficient carrier gate fails at `ell=23203` and passes at `ell=23204`. This is a reproducible boundary of the implemented sufficient criterion, not a manuscript parameter choice or a claimed optimal threshold.

A successful `LargeBandPhaseScaleCertificate` exposes only the consequent scalar bounds

- `phaseError <= 4/S_*`, and
- `phaseConstant(M) * phaseError <= phaseConstant(M) * 4/S_*`.

This closes a representation gap between the asymptotic Lean theorem and the robust binary64 `DyadicChart` path without changing `charts.py` or weakening its numerical safeguards.

## Independent verification

`tests/test_phase_large_band_scale.py`:

- independently checks the `h=1/200` carrier boundary using the log-space inequality `4 log2(ell) <= ell/400`, rather than recomputing the production integer powers;
- verifies the certificate still works when binary64 `2^-ell` has underflowed to zero;
- checks the simplified normal envelope directly against `phase_constant(M) * 4/ell^2`;
- exercises invalid band, `h`, and `M` inputs;
- asserts all paper-exact truth flags remain false.

## What this does not prove

This module does **not** provide the paper-exact Proposition 5.5 / Section 5 background field, does not certify the `LocalBaseBounds` C1/C2 hypotheses on any actual slow box, does not prove representative shear orthogonality or actual/reference normal closeness, and therefore does not by itself establish Eqs. (7.9)–(7.11) uniformly on the active family. It also does not select stress-cone amplitudes or construct curl-realized waves.

Accordingly:

- `actual_base_fields_verified = false`;
- `uniform_eq_7_9_to_7_11_verified = false`;
- `paper_exact_velocity_available = false`.
