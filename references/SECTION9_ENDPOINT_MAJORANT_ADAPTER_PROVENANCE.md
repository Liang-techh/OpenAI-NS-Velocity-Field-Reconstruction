# Section 9 Eq. (9.18) -> Section 10 endpoint-majorant adapter provenance

Status: **formal-structure**. This increment derives one endpoint-majorant
coefficient from theorem-facing Section 9 data and the already-fixed Section 10
late-time geometry. It does **not** construct the manuscript correction
sequence, prove Eq. (9.18), identify the actual Eq. (9.21) residual provider, or
prove the endpoint limits/Borel extension.

## Gap closed

Before this adapter, `section9_endpoint_bridge.py` correctly required a
separately justified uniform physical-time estimate

`sup_x |partial_t D^n R(t,x)| <= C (1-t)^(-alpha)`

but `C` and `alpha` were still direct theorem inputs. That fail-closed API is
useful, but it leaves a free conversion step between the Section 9 Eq. (9.18)
residual envelope and the Section 10 endpoint theorem.

`section9_endpoint_majorant_adapter.py` removes that conversion freedom for the
bounded case. A caller supplies a genuinely uniform Eq. (9.18) envelope for
derivative order `m=n+1`,

`|R(u^[j],p^[j])|_m <= C q^beta (1+|log q|)^P + E`,

with theorem provenance for both `C` and a uniform flat-remainder bound `E`.
The adapter computes `beta = h*sigma_j-K_m` from the existing exact Section 9
exponent ledger and derives the endpoint coefficient analytically.

## Official physical-domain bound

This step reuses, rather than replaces, the landed Section 10 geometry. On the
official late plateau `t>=3/4` and the fixed spatial support `|z|<=1/4`, the
Eq. (4.1) argument already encoded in `section9_physical_window.py` gives

`0 < q <= (1-t)+z^2 <= 1/4+1/16 = 5/16`.

The witness must certify that its uniform Eq. (9.18) theorem is valid on a
strict q-domain extending beyond `5/16`. No arbitrary time window is introduced.

## Analytic power-log maximization

For `beta>0`, put `s=-log(q)`. On `0<q<=5/16`,

`q^beta (1+|log q|)^P = exp(-beta s)(1+s)^P`

with `s>=-log(5/16)`. The logarithmic derivative is

`-beta + P/(1+s)`,

so the maximum is attained either at the left endpoint or at
`s=P/beta-1`. The implementation evaluates exactly that stationary-point rule.
It also admits the constant case `beta=0, P=0`. It rejects `beta<0` and
`beta=0, P>0`, because those cases do not yield a bounded endpoint majorant by
this argument.

The resulting theorem-facing endpoint witness has `alpha=0`, hence an
integrable physical-time majorant on `[3/4,1)`. It is then passed through the
existing `admit_section9_uniform_endpoint_witness` and official Section 10
localization-transfer gate.

## Fail-closed source semantics

`Section9UniformResidualEnvelopeWitness` requires:

- endpoint degree `n`, Section 9 stage, and spatial-window identity;
- the existing exact `h*sigma_j-K_m` exponent inputs;
- theorem/certified evidence for the leading constant and a **uniform**
  flat-remainder bound, from the same evidence class;
- a strict q-validity upper endpoint greater than `5/16`;
- a stable `(source_id, source_revision)`;
- explicit theorem assertions that the Eq. (9.18) envelope is uniform, the flat
  remainder bound is uniform, and the theorem covers the official Section 10
  support.

These booleans are theorem inputs, not facts inferred by Python. Sampled grids
cannot promote themselves through this interface merely by having small values.

## Regression

`tests/test_section9_endpoint_majorant_adapter.py` checks:

- exact exponent bookkeeping and the derived `alpha=0` bridge;
- the fixed `5/16` Section 10 q-bound;
- the analytic supremum against an independent dense-grid regression;
- the admissible constant case `beta=0,P=0`;
- rejection of negative/non-decaying power-log envelopes;
- rejection when the theorem q-domain does not cover `5/16`;
- rejection of mixed evidence classes or missing theorem assertions;
- preservation of all fail-closed truth flags.

The dense grid is only an independent regression check. The production
coefficient is obtained from the analytic stationary-point calculation, not from
sampling or fitting.

## Remaining Issue #4 boundary

This adapter is reusable analytic infrastructure for the true `t=1` path, but
it does not make the source real. Still required are:

1. materialize the genuine Section 9 correction sequence and Eq. (9.21) residual
   source on the late-time support;
2. obtain source-linked uniform Eq. (9.18) constants/losses/flat-remainder
   bounds for every required spacetime derivative order;
3. bind the resulting machine-derived endpoint-majorant ladder to that exact
   residual revision and construct the actual locally uniform endpoint limits;
4. complete the infinite Borel right-extension/all-order smoothness argument;
5. only then reconstruct the compact smooth forcing and independently check
   closure/smoothness, divergence-free structure, support, finite energy, and
   the blow-up path.

Accordingly `actual_section9_sequence_verified=false`,
`source_majorant_derived_from_actual_residual_verified=false`,
`endpoint_limit_constructed=false`, and
`paper_exact_velocity_available=false` remain mandatory.
