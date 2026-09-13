# Hierarchy-owned positive-axis forcing second-parameter provenance

## Scope

This increment adds one analytic derivative-composition layer for Section 5. The exact Eq. (5.7) forcing value and its first eta derivative were already landed. The strict-lower repaired hierarchy now also owns the second eta derivative of `actualLowerSource`. The new bridge composes those two existing layers to expose

`(f_n, partial_eta f_n, partial_eta^2 f_n)`

for the same strong repaired lower-history hierarchy.

## Paper/formula binding

The forcing is the displayed six-vector from the positive-axis first-order system (Eq. (5.7)):

- rows 1--3 are zero;
- row 4 is `2 xi * pressureSource`;
- row 5 is `2 * actualLowerSource.angular`;
- row 6 is `2 * actualLowerSource.axial - 4 X (eta/ell) pressureSource`, with `X=xi^2` and `ell=1-2 h eta^2`.

`pressureSource = C^-2 * pressureProduct - (Omega/X)/2` remains the already-landed exact linear combination, so its first and second eta derivatives are obtained by applying the same linear map to the hierarchy-owned source jets.

Only the explicit scalar geometry is differentiated in this increment. For `g(eta)=eta/ell(h,eta)`, production uses

`g' = 1/ell + 4 h eta^2/ell^2`

and

`g'' = 12 h eta/ell^2 + 32 h^2 eta^3/ell^3`.

Thus the new sixth-row second derivative is assembled analytically as

`2 source_axial'' - 4 X (g'' p + 2 g' p' + g p'')`.

## Ownership and fail-closed boundary

The public bridge accepts only `Section5LowerHistoryPhiFourthMixedHierarchy`. It does not accept an independent forcing derivative, source derivative, beta/Omega table, generic cutoff, or fitted coefficient family. The hierarchy must therefore own the repaired fourth-mixed phi / fifth-mixed U data already required by the second-eta strict-lower-source bridge. Missing genuine order-zero strong data continue to fail closed, including at `xi=0`; this increment does not synthesize Issue #1 data.

The forcing value and first derivative delegate to the previously landed hierarchy-owned forcing bridge and remain authoritative there. Only the second derivative is new.

## Regression

`tests/test_background_repaired_history_forcing_second_parameter.py` uses the actual Lemma 5.2 compact five-bump repair for the positive order in a recursive order-2 history. It checks points on the axis and inside both compact-repair regions. The new value and first derivative must equal the previously landed bridge exactly, while `partial_eta^2 f_n` is cross-checked against a centered eta derivative of that prior analytic first-derivative path. Finite differences are test-only oracles and are not used in production.

The regression also verifies that a hierarchy lacking leading fourth-mixed phi ownership fails closed at the axis.

## Truth boundary

Status remains `formal-structure`. `paper_exact_velocity_available=false` and `full_reconstruction=false`.

This increment does **not** establish `partial_eta W_n^(1)`, Picard convergence, final recursive coefficient materialization, the paper's recursive cutoff-scale schedule, Proposition 5.3 all-jets residual decay, or paper-exact velocity. Caller/test leading fourth/fifth-mixed data are not paper-exact substitutes for the genuine Issue #1 leading profile.
