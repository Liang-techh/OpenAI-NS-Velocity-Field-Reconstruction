# Section 5 hierarchy-owned axial preceding-diffusion third-eta provenance

## Scope

The Section 5 strict-lower PositiveAxis source contains an axial preceding-
diffusion contribution

`Z_(b-D) (Z_b U_(n-1))`

with the axial paper power `-1/2-h`.  PR #280 landed the exact analytic
third-eta preceding-diffusion primitive and the hierarchy-owned angular bridge.
This increment adds only the corresponding fail-closed axial bridge.

The implementation baseline is repository `main`
`aa3f0ee90d871fbbcfb5820471aa0f3dd7b7bb7e`, which already includes PR #280.

## Hierarchy ownership

`hierarchy_owned_axial_preceding_diffusion_third_parameter_jet` accepts only a
`Section5LowerHistoryFifthMixedHierarchy`.  For positive coefficient order `n`
it queries the hierarchy's authoritative repaired fifth-mixed `U_(n-1)` jet,
re-types the identical mixed-jet rows into the generic
`ProfileFifthMixedJet` carrier, and calls the landed exact analytic
`preceding_diffusion_third_parameter_jet` with power `-1/2-h`.

The fifth-mixed U provider is the already-landed Lemma 5.2 compact five-bump
repair layer.  Its fourth- and third-mixed projections are checked by the
hierarchy for exact coherence.  The new bridge accepts no caller-maintained
source derivative table, no separately supplied preceding-diffusion jet, no
sampled fit, and no generic cutoff.

Genuine order-zero strong profile data and the unrepaired fifth-mixed provider
remain upstream contracts.  This bridge therefore reduces a solver ownership
seam but does not promote those contracts to paper-exact data.

## Verification

`tests/test_background_repaired_history_axial_preceding_diffusion_third_parameter.py`
uses an actual `Lemma52MomentRepair.from_intervals(...)` five-bump repair and
evaluates inside active U-repair support at `X = 0.5 * 2.17^2`.

The regression checks that:

1. the hierarchy-owned repaired `U_1` fifth-mixed jet differs from the
   unrepaired base fixture inside the active bump;
2. the returned value/first/second eta rows agree exactly with the already
   landed second-eta preceding-diffusion implementation;
3. the new third eta row agrees with centered eta differentiation of that older
   analytic second-eta path, with finite differencing used only as a test
   oracle; and
4. missing strict-lower fifth-mixed history and wrong hierarchy types fail
   closed.

The polynomial moment/base providers in the regression are test fixtures only
and are not paper coefficient data.

## Truth boundary

Status remains Stage-2 `formal-structure` with `full_reconstruction=false` and
`paper_exact_velocity_available=false`.

After this increment both angular and axial third-eta preceding-diffusion rows
have hierarchy-owned bridges.  This still does not establish the complete
hierarchy-owned `partial_eta^3 actualLowerSource`: in particular the third-eta
`Omega_(n-1)/X` chain requires a stronger regular-flux/Omega derivative layer.
Consequently hierarchy-owned `partial_eta^3 f_n`, `partial_eta^3 W_n^(0)`,
`partial_eta^2 W_n^(1)`, Picard convergence/final coefficient materialization,
the paper recursive cutoff-scale schedule, uniform `C[j,m]` bounds,
Proposition 5.3 all-order residual decay, full reconstruction, and paper-exact
velocity remain open.

Uniform bounds, recursive cutoff selection, and all-order convergence are
explicitly outside this increment and remain for the separate convergence
workstream.  Green CI for this bounded bridge must not be interpreted as an
all-order or paper-exact result.
