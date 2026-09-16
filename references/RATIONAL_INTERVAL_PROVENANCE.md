# Shared exact-rational interval arithmetic

Status: implementation utility only. This refactor does not add a mathematical constructor, certify the full pressure integral, or change `paper_exact_velocity_available=false` / `full_reconstruction=false`.

## Scope

NS003 consolidates exact-rational interval primitives used by the certified outgoing sigma, tail-debt, and release-lag enclosures. The previous implementations duplicated positive-control validation and dyadic grid helpers, while `outgoing_release_lag_enclosure.py` imported several private helpers from `outgoing_tail_debt_enclosure.py`.

The new `rational_interval.py` owns:

- ordered `Fraction` interval validation;
- exact interval addition, subtraction, sign-safe multiplication, and division by a proved strictly positive denominator interval;
- endpoint propagation through exact rational-valued monotone functions;
- outward dyadic floor/ceil rounding with optional proved clipping bounds;
- common positive rational/cap validation and dyadic-step selection.

`outgoing_sigma_enclosure.py` continues to re-export `RationalInterval`, preserving its existing public import surface. Its positive exponential helper is now a public `validated_exp_positive_on_0_3` entry point. Tail-debt and release-lag both use that same implementation, so release-lag no longer reaches into tail-debt for private grid/exponential helpers.

## Mathematical semantics

All operations use exact `fractions.Fraction` endpoints. Multiplication considers all four endpoint products, so sign changes are not treated as a positive-only special case. Positive interval division explicitly rejects denominators with lower endpoint `<= 0`; near-zero positive denominators are enlarged by exact reciprocal endpoint reversal rather than regularized. Monotone propagation evaluates only exact endpoints and rejects a direction inconsistent with the returned order. Dyadic rounding always floors the lower endpoint and ceils the upper endpoint before any independently justified clipping bounds are intersected.

The outgoing sigma, tail-debt, and release-lag formulas and tolerance stopping criteria are otherwise unchanged. The tail ratio `rho/(1-rho)` and positive product in the debt path are routed through the shared arithmetic. No old floating evaluator becomes certification evidence.

## Verification

Independent utility regression covers sign-changing products, exact zero, a denominator interval as small as `1/10^12`, rejection of a denominator touching zero, increasing/decreasing monotone endpoint propagation, outward dyadic rounding, and bool/zero/negative cap rejection.

A local isolated execution of the same utility tests before repository submission reported `5 passed in 0.06s`. The execution container could not resolve `github.com`, so a clean repository checkout was not available there. Repository GitHub Actions on the PR head are authoritative for the existing outgoing sigma/tail-debt/release-lag regressions and broader suites. No CI result is claimed until those runs complete.

## Boundary

This change is an arithmetic refactor and shared safety layer. It does not certify decayHold, transition geometry, clockWeight, pressure jets, Stage 1 fixed-point closure, or any later reconstruction stage. Passing these tests cannot promote paper-exact or full-reconstruction status.
