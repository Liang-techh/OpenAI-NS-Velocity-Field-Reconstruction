# Section 5 hierarchy-owned lower-source second-eta provenance

Status: **Stage 2 / formal-structure only**.

## Scope

PR #240 composes already-landed Section 5 derivative layers into the first
hierarchy-owned analytic second eta jet of the complete strict-lower
`actualLowerSource` used by the PositiveAxis Eq. (5.7) forcing.

The new bridge is
`hierarchy_owned_lower_source_second_parameter_jet(...)` in
`background_repaired_history_source_second_parameter.py`.

It does not introduce a new paper formula. It reuses:

- the landed value/first-eta `actualLowerSource` bridge;
- hierarchy-owned repaired fourth-mixed `phi` data from Lemma 5.2;
- hierarchy-owned repaired fifth-mixed `U` data and its lower projections;
- fourth-mixed `beta=V/X` derived only through the analytic Eq. (5.2) mapper;
- the analytic second-eta preceding-diffusion operator; and
- the hierarchy-owned analytic second-eta Eq. (5.6) `Omega/X` bridge.

## Paper / implementation map

- Lemma 5.2 / Appendix A: compact moment repair supplying the strong positive-order history.
- Eq. (5.2): derives the regular radial flux `beta_n=V_n/X`, including the already-landed positive-order `D+lambda_n` correction.
- Eqs. (5.3)-(5.6): define the strict-lower angular, axial, pressure-product, preceding-diffusion, and regular `Omega/X` source rows.
- Eq. (5.7): consumes this source downstream in the PositiveAxis first-order system.

The new implementation differentiates the existing strict-lower convolution
rows once more in eta by the ordinary product rule, subtracts the already-landed
second-eta preceding-diffusion rows, and uses the already-landed hierarchy-owned
`partial_eta^2(Omega/X)` value. The pre-existing first-parameter bridge remains
authoritative for both the source value and its first eta derivative.

## Verification

Regression uses the actual `Lemma52MomentRepair.from_intervals(...)` compact
five-bump geometry and the strong repaired order-1 coefficient inside a
recursive order-2 history.

At `X=0` and points inside both compact-repair regions:

1. the new source value is required to equal the previously landed source value exactly;
2. the new first eta derivative is required to equal the previously landed analytic first-parameter bridge exactly; and
3. each component of the new second eta derivative is cross-checked against a centered eta derivative of that older analytic first-derivative path.

Finite difference is test-only. Production contains no numerical eta
differentiation.

A hierarchy that owns the preceding U-fifth layer but not the stronger phi-fourth
layer fails closed, including at the axis. This prevents an axis-zero shortcut
from hiding missing Issue-#1 strong data.

## Truth boundary

The genuine order-zero fourth/fifth mixed leading profile remains an upstream
Issue #1 input. Positive-order unrepaired base jets and the moment/patch eta jets
also remain inputs to the compact repair.

Therefore this increment does **not** claim any of the following:

- paper-exact leading or positive-order coefficient profiles;
- hierarchy-owned `partial_eta^2 f_n`;
- hierarchy-owned `partial_eta W_n^(1)` or a subsequent Picard iterate;
- Picard convergence or final recursive coefficient materialization;
- the paper's recursively chosen cutoff-scale schedule;
- Proposition 5.3 all-jets truncation/residual decay; or
- paper-exact reconstructed velocity.

`paper_exact_velocity_available=false` and `full_reconstruction=false` remain
mandatory. Caller-supplied strong leading jets, test polynomial families,
generic cutoffs, sampled derivatives, and finite-order numerical fits are not
promoted to paper-exact data by this bridge or its regression tests.
