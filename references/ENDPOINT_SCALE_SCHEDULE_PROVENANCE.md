# Section 10 endpoint Spatial-Borel scale-schedule provenance

Status: **formal-structure / constructive admissible scale witness**. This file does not promote Stage 7 or the final velocity/force to paper-exact.

## Pinned source

Official repository: `openai/NavierStokesAndEuler`

Pinned commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Primary file: `NavierStokes/SpatialBorelExtension.lean`.

Relevant definitions/theorems:

- `template m j a = localizedTerm m 1 j a`;
- `exists_template_bound`, `templateBound`, and `template_deriv_le` obtain a global nonnegative derivative bound for each compactly supported smooth template;
- `boundSum a ha j = sum_{m<j} sum_{k<j} templateBound (ha j) m j k`;
- `localScale a ha j = Classical.choose (exists_nat_gt (2^j * boundSum a ha j))`;
- `scale = DiagonalScale.doublingEnvelope localScale`;
- `boundSum_le_scale` gives `boundSum(j) <= 2^{-j} scale(j)`;
- `localized_derivative_tail_bound` combines the scaled-template derivative estimate with `m<j`, `k<j` to obtain the geometric tail `<=2^{-j}`;
- `majorant` is exactly this geometric `2^{-j}` quantity once `j>max(m,k)`, with only finitely many lower-degree exceptions;
- `majorant_summable` and `localized_derivatives_summable` use that eventual geometric majorant to prove summability of every fixed localized derivative series;
- `derivative_tail_bound_on_compact`, `localizedExtension_contDiff`, and `extension_contDiff` then use these bounds to prove all-order local smoothness of the Borel series.

The downstream consumer is `NavierStokes/SpacetimeGluing.lean`; `NavierStokes/CandidateFromLimits.lean` supplies endpoint jets only after assuming locally uniform limits of every full spacetime derivative of the actual Navier-Stokes residual.

## Executable increment

`src/openai_ns_reconstruction/endpoint_scale_schedule.py` makes the scale and eventual-tail arithmetic above constructive:

1. a caller supplies `template_bound(m,j,k)`, intended to be a valid global bound for the official localized template derivative;
2. `SpatialBorelScaleSchedule.bound_sum(j)` forms the exact finite triangular sum over `m,k<j`;
3. instead of Lean's opaque `Classical.choose`, `local_scale(j)` deterministically chooses the least natural strictly greater than `2^j * bound_sum(j)`;
4. the existing `DoublingEnvelope` supplies the official doubling recurrence;
5. `DegreeScaleCertificate` checks the exact rearranged inequality `bound_sum(j) <= 2^{-j} scale(j)`;
6. `TailScaleCertificate` checks the exact supplied-bound version of `(scale^k/scale^j) templateBound <= 2^{-j}` for each `m,k<j`;
7. `DerivativeTailSumCertificate` records the infinite omitted-series budget: if degree `N` is retained and `N>=max(m,k)`, every omitted `j>N` lies in the geometric regime, hence
   `sum_{j>N} ||D^k localizedTerm_j|| <= sum_{j>N} 2^{-j} = 2^{-N}`;
8. `borel_extension_from_template_bounds` connects this schedule to the already-landed full-spacetime-to-normal jet adapter and Taylor-Borel evaluator.

Scale and tail inequalities are evaluated with `fractions.Fraction`: integer/Fraction inputs remain exact, and finite float-like inputs are converted to their exact binary64 rational value before schedule arithmetic. The new all-omitted-degrees budget is therefore represented as the exact rational `1/2^N`, not inferred from sampled convergence or a fitted decay rate.

## Independent checks

`tests/test_endpoint_scale_schedule.py` independently recomputes the finite bound sums, the least strict-natural witness, the doubling recurrence and the localized derivative ratio. For the summed derivative tail, it checks a finite block using the direct scaled-template bounds, closes the remaining infinite part by an independently evaluated geometric remainder, and separately verifies the closed form `sum_{j>N}2^{-j}=2^{-N}`. It also checks that the scale path used by the existing Borel evaluator agrees with the schedule and exercises negative/nonfinite/malformed inputs as fail-closed cases.

## Deliberate boundary

This increment **does not construct `templateBound` from the actual residual jets**. In Lean those bounds follow from smoothness plus compact support of the localized template. The executable repository still lacks the actual completed local velocity/pressure from upstream stages, the full derivative family of their activated Navier-Stokes residual, locally uniform `t -> 1-` limits at every order, and analytic all-order template derivative bounds for those true endpoint jets.

Consequently, the new `2^{-N}` tail certificate is an exact conditional majorant statement, not by itself a proof that the repository's future force is `C^infinity` through `t=1`. A caller-supplied bound provider is not a paper-exact certificate by itself, finite sampled derivative maxima are not accepted as a substitute, and the resulting extension must not be described as a proved smooth global force. `paper_exact_velocity_available` remains `false` and Stage 7 remains `formal-structure`.
