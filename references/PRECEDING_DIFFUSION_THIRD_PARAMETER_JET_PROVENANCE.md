# Section 5 preceding-diffusion third parameter-jet provenance

## Scope

The PositiveAxis strict-lower source contains

`precedingDiffusion = Z_(b-D) (Z_b F_(n-1))`.

The repository already owns the value, first eta derivative, and second eta
derivative of this term.  This increment adds only the exact analytic third eta
derivative and a narrow angular bridge from the strong repaired Section-5
hierarchy.

The implementation baseline is repository `main`
`67ab5ecd32f80225df8bce2760b943c20a7d3d7b`, which includes PR #278 and its
hierarchy-owned repaired fifth-mixed `phi_n` provider.

## Paper-derived derivative structure

Write the inner operator as

`M = Z_b F = N / ell`,

where `ell = 1 - 2 h eta^2`.  The outer operator is

`precedingDiffusion = Z_(b-D) M = Q / ell`.

To compute `partial_eta^3(Q/ell)` exactly, the inner row requires
`M^(0)..M^(4)` and `partial_X M^(0)..partial_X M^(3)`.  Relative to the landed
fourth-mixed profile jet, this introduces exactly three total-order-five input
fields:

- `F_XXetaetaeta`,
- `F_Xetaetaetaeta`, and
- `F_etaetaetaetaeta`.

`background_preceding_diffusion_third_parameter_jet.py` consumes the landed
`ProfileFifthMixedJet`, delegates the value/first/second rows to
`preceding_diffusion_second_parameter_jet`, and evaluates only the new third
row using explicit analytic derivatives of `N`, `partial_X N`, and
`1 / ell`.  Production performs no finite differencing, sampled derivative
fit, or generic-cutoff substitution.

## Hierarchy-owned angular bridge

`hierarchy_owned_angular_preceding_diffusion_third_parameter_jet` accepts only
`Section5LowerHistoryPhiFifthMixedHierarchy`.  For positive coefficient order
`n`, it queries the authoritative hierarchy-owned repaired
`phi_(n-1)` fifth-mixed jet before evaluating the angular paper power
`-1-h`.  The query therefore fails closed when the genuine strong strict-lower
history is absent, including at the axis; no caller-maintained derivative table
can bypass that ownership check.

The positive-order repaired source used by this hierarchy is the actual
Lemma 5.2 compact five-bump adapter.  Its unrepaired base fifth-mixed profile,
normalization `C`, and fifth-order moment/patch eta jets remain explicit
upstream contracts.

## Verification

`tests/test_background_preceding_diffusion_third_parameter_jet.py` verifies:

1. exact identity of the returned value/first/second rows with the already
   landed second-parameter implementation;
2. the new analytic third derivative against a centered eta derivative of the
   older analytic second-derivative path, with finite differencing used only as
   a test oracle;
3. sensitivity to each of the three newly required total-order-five fields;
4. fail-closed behavior for weaker and non-finite stronger inputs; and
5. composition through a real `Lemma52MomentRepair.from_intervals(...)`
   five-bump repair inside active E-repair support, where the hierarchy-owned
   repaired fifth-mixed `phi_1` is passed into the new angular operator.

Polynomial profiles in the derivative regression are test fixtures only and
are not promoted to paper coefficient data.

## Truth boundary

Status remains Stage-2 `formal-structure` with
`full_reconstruction=false` and `paper_exact_velocity_available=false`.

This increment establishes the hierarchy-owned **angular** third-eta
preceding-diffusion row only.  It does not establish the corresponding axial
third-eta row, the full hierarchy-owned `partial_eta^3 actualLowerSource`,
`partial_eta^3 f_n`, `partial_eta^3 W_n^(0)`, `partial_eta^2 W_n^(1)`, Picard
convergence, finalized positive-order coefficient profiles, the paper's
recursive cutoff-scale schedule, Proposition 5.3 arbitrary-order residual
decay, full reconstruction, or paper-exact velocity.

Issue #1 remains the source of genuine leading strong profile data.  Green CI
for this bounded analytic bridge must not be interpreted as evidence that those
upstream data or the final velocity field have been reconstructed.
