# Positive-axis forcing third-parameter primitive provenance

## Scope

This artifact is a bounded Stage-2 / Issue #2 solver prerequisite for Section 5.  It differentiates only the displayed Eq. (5.7) positive-axis forcing one additional time in the similarity parameter `eta`.

The current repaired lower-history hierarchy on `main` owns `f_n`, `partial_eta f_n`, and `partial_eta^2 f_n`, but the exact second-eta Picard recursion primitive landed in PR #272 exposes a genuine next requirement: the previous iterate needs a third eta derivative, and for `W_n^(0)=G f_n` that in turn requires `partial_eta^3 f_n`.

## Exact identity implemented

For the sixth forcing row write

`g(eta) = eta / ell`, with `ell = 1 - 2 h eta^2`,

and `p = pressureSource`.  The forcing contribution is `-4 X g p`.  The implementation uses

- `g' = 1/ell + 4 h eta^2 / ell^2`,
- `g'' = 12 h eta / ell^2 + 32 h^2 eta^3 / ell^3`,
- `g''' = 12 h / ell^2 + 192 h^2 eta^2 / ell^3 + 384 h^3 eta^4 / ell^4`,

so

`partial_eta^3(g p) = g''' p + 3 g'' p' + 3 g' p'' + g p'''`.

Rows four and five are differentiated directly from `2 xi pressureSource` and `2 actualLowerSource.angular`; the first three rows remain zero.

## Input boundary

`PositiveAxisSourceThirdParameterJet` is an explicit analytic source-jet contract.  It is **not** claimed to be hierarchy-owned or paper-exact.  The current main branch still lacks a strong repaired lower-history construction of `partial_eta^3 actualLowerSource`; that missing ownership remains the next Section-5 derivative debt.

The production path contains no eta finite differences, sampled derivative tables, numerical fit, generic cutoff substitution, or caller-supplied recursive coefficient family labelled as exact.  Centered differences appear only in regression as an independent oracle obtained by differentiating the already-displayed second-eta forcing formula.

## Verification

`tests/test_background_positive_axis_forcing_third_parameter.py` checks:

1. the analytic third derivative against a centered derivative of an independently assembled exact second-eta forcing formula for cubic source data;
2. the isolated effect of a pressure-source third derivative on the two rows where it is permitted to enter;
3. that the axis does not erase nonpressure third-source data; and
4. fail-closed behavior for malformed source jets and invalid coordinates.

## Truth boundary

Status remains `formal-structure`.

This artifact does **not** establish hierarchy-owned `partial_eta^3 actualLowerSource`, hierarchy-owned `partial_eta^3 f_n`, hierarchy-owned `partial_eta^2 W_n^(1)`, Picard convergence, finalized positive-order coefficients, the paper recursive cutoff-scale schedule, Proposition 5.3 residual decay, or paper-exact velocity.  `full_reconstruction=false` and `paper_exact_velocity_available=false` remain mandatory.
