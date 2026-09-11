# Lemma 5.2 functional repaired-profile provenance

Status: **formal-structure only**.

This increment does not materialize the profile-dependent Section-5 hierarchy and
does not change `paper_exact_velocity_available=false`.

## Pinned source

Official Lean source commit:
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Paper locations: Lemma 5.2, Eqs. (5.14)-(5.16), with downstream use in
Eq. (5.15).

The repository already had two lower layers:

- `background_moment_repair.py` realizes the two fixed `U_n` bumps and three
  fixed `E_n` bumps and solves the two constant moment systems at one `eta`;
- `background_moment_repair_jets.py` propagates supplied ordinary eta
  derivatives of the five unrepaired moments and of `p(eta)=e_* f(eta)` through
  the exact quotient/Leibniz recurrences, producing derivatives of the repair
  coefficients.

## What this increment adds

`background_moment_repair_profile.py` turns those repair coefficient jets into
actual eta-dependent repaired coefficient functions:

`U_n(X,eta) = U_tilde_n(X,eta) + sum_j alpha_j(eta) b^U_j(sqrt(2X))`,

`E_n(X,eta) = E_tilde_n(X,eta) + sum_j beta_j(eta) b^E_j(sqrt(2X))`.

Because the radial bumps are fixed and independent of `eta`, the adapter also
implements

`d_eta U_n = d_eta U_tilde_n + sum_j alpha'_j b^U_j`,

which is precisely the derivative required by the landed Eq. (5.2) radial-flux
formula and by `background_extension.reconstruct_eq_5_15`.

The adapter can expose the repaired functions through the repository
`LeadingProfile` API, but `paper_exact` is hard-coded to `False`; a pressure
profile is never inherited silently from the unrepaired input.

## Independent regression checks

The new tests do not merely inspect the coefficient arrays returned by the same
linear solve.  For a nonconstant eta-dependent moment/patch datum they:

1. integrate the actual function-level differences `U_repaired-U_base` and
   `E_repaired-E_base` with independent SciPy adaptive quadrature and verify all
   five Eq. (5.16) corrected moments vanish to numerical tolerance;
2. compare the analytic repaired `d_eta U_n` against an independent centered
   parameter finite difference of the repaired `U_n` function;
3. verify exact equality with the base profile outside the five compact radial
   supports; and
4. verify that missing eta derivatives or a zero center patch factor fail
   closed.

These are executable cross-checks of the implication from supplied jets to a
functional compact repair.  They are not a substitute for the manuscript's
uniform analytic estimates.

## Truth boundary / remaining blocker

The adapter still accepts the unrepaired moment jets and `e_* f` jets from its
caller.  It does **not**:

- derive those jets from the actual converged Eq. (5.7) positive-order
  coefficient;
- prove the manuscript's uniform nonvanishing lower bound for `e_* f` on the
  eta interval;
- certify analytic-strip bounds or all eta derivatives from finite samples;
- construct the true recursive `A0/A1/f_n` data blocked on Issue #1; or
- prove next-order support/stress closure or Proposition 5.3 residual flatness.

The next profile-dependent step remains to materialize the genuine Eq. (5.7)
coefficient from the Section-4 leading profile, derive its five unrepaired
moment functions/eta jets, and feed those certified data into this adapter
before invoking the existing Eq. (5.15) reconstruction and next-order
recursion.
