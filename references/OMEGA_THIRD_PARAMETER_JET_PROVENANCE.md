# Section 5 hierarchy-owned third-eta `Omega/X` provenance

## Scope

This increment closes the next derivative seam on the route from the repaired strict-lower coefficient history to a hierarchy-owned third eta forcing jet. The implementation baseline is repository `main` commit `1ea24db9524d34a191ebf07d3e22a7b5a1ecdd6f`, which already includes PR #289 and therefore owns compact Lemma 5.2 repaired sixth-mixed `U_n` jets and analytically derived fifth-mixed `beta_n=V_n/X` jets.

The new primitive evaluates value through `partial_eta^3(Omega_k/X)` for Eq. (5.6). The landed second-eta implementation remains authoritative for value, first eta derivative, and second eta derivative; only the new third derivative is assembled here.

## Exact analytic differentiation

All five Eq. (5.6) rows are differentiated analytically:

- the regular `T` row of `beta_k`;
- the quadratic regular-flux convolution;
- the `U_i Z(beta_j)` convolution;
- the radial viscosity `-2(2 beta_X + X beta_XX)`;
- the shifted nested viscosity `-Z_(b-D) Z_b(X beta_{k-1}) / X`.

Every quotient by `L=1-2h eta^2` is differentiated through the exact identity `L q = N`. For derivative order `m`, because `L` is quadratic,

`q^(m) = [N^(m) - m L' q^(m-1) - binom(m,2) L'' q^(m-2)] / L`.

Thus production does not use finite differences, fitted coefficients, sampled `V/X` division, or a generic cutoff.

The nested shifted-viscosity row needs the inner quotient through eta order four and its X derivative through eta order three. The fifth-mixed regular-flux fields landed in PR #289 — `beta_XXetaetaeta`, `beta_Xetaetaetaeta`, and `beta_etaetaetaetaeta` — are exactly sufficient. No seventh-mixed U layer is introduced.

## Hierarchy ownership

`hierarchy_owned_omega_third_parameter_jet` accepts only `Section5LowerHistorySixthMixedHierarchy`. For every required lower order it obtains `RegularFluxFifthMixedJet` from `beta_fifth_mixed_jet`, which derives beta only through analytic Eq. (5.2), and obtains the axial third-eta row as an exact projection of the same hierarchy-owned sixth-mixed U jet.

Missing strong history, including the genuine Issue #1 leading coefficient when it is not supplied, fails closed. No caller-maintained beta table or hand-filled `SourceJet` is accepted by this bridge.

## Verification

`tests/test_background_omega_third_parameter_jet.py` uses analytic polynomial beta/U fixtures. At positive order one it:

- verifies value, first eta derivative, and second eta derivative exactly equal the landed second-eta implementation;
- exercises the nested shifted-viscosity row;
- cross-checks the new third eta derivative against centered differentiation of the landed analytic second-eta row, with finite differences confined to the test oracle;
- checks regular evaluation on the axis `X=0`.

The polynomial data are regression fixtures only and are not paper coefficients.

## Truth boundary

Status remains Stage-2 `formal-structure`, with `full_reconstruction=false` and `paper_exact_velocity_available=false`.

This increment does not establish the complete `partial_eta^3 actualLowerSource`, hierarchy-owned `partial_eta^3 f_n`, `partial_eta^3 W_n^(0)`, `partial_eta^2 W_n^(1)`, Picard convergence, final coefficient materialization, uniform `C[j,m]` bounds, recursive cutoff convergence, Proposition 5.3 all-order residual decay, full reconstruction, or paper-exact velocity.

Uniform bounds, recursive cutoff scheduling, and all-order convergence remain with Agent 7's separate convergence workstream.
