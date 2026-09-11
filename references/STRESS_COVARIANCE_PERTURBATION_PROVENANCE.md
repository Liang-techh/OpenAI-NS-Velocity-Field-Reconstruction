# Section 7 actual-vs-reference covariance perturbation provenance

Status: **formal-structure**. This increment does not upgrade Stage 3/4 or the
repository-wide paper-exact gate.

## Published construction mapped here

Paper: OpenAI, *Finite Time Blowup for Navier-Stokes* (September 2026),
Section 7, Proposition 7.5, especially Eqs. (7.27)-(7.29), printed p. 83 in
the inspected PDF.

The paper writes the actual pulse covariance columns in the frozen `(N,K)`
frame as

`H_sigma = h_sigma (-A_c N - sigma u_* K + e_sigma)`

with a uniform error `|e_sigma| <= C S_*^(-1/2)` in Eq. (7.28). It then uses
the strict order-zero cone margin and smallness of these errors to obtain the
nondegeneracy, inverse and positive-coefficient estimates in Eq. (7.29).

The landed `stress_cone.py` already implements only the zero-error signed
reference problem. This increment adds the finite-dimensional implication from
a *previously certified* column-error bound to persistence of invertibility and
positivity. It does not construct the true `H_sigma` integrals or prove the
paper's analytic `C S_*^(-1/2)` estimate.

## Executable perturbation estimate

Let the normalized reference columns be

`v_- = (-a,-b)`, `v_+ = (-a,+b)`,

with positive column scales `s_-`, `s_+`, and let the supplied normalized errors
satisfy `|e_-|, |e_+| <= delta`. The actual matrix represented by the module is

`H = [ s_- (v_- + e_-),  s_+ (v_+ + e_+) ]`.

Writing `r=sqrt(a^2+b^2)`, determinant multilinearity gives

`|det([v_-+e_-,v_++e_+]) - det([v_-,v_+])|`
`<= 2 r delta + delta^2`.

Since the reference determinant is `-2ab`, the strict scalar condition

`2ab - (2 r delta + delta^2) > 0`

preserves the determinant sign and yields the explicit lower bound

`|det H| >= s_- s_+ [2ab - (2 r delta + delta^2)]`.

For the spectral inverse norm, the 2x2 singular-value identity gives

`||H^-1||_2 = sigma_max(H)/|det H| <= ||H||_F/|det H|`,

while

`||H||_F <= (r+delta) sqrt(s_-^2+s_+^2)`.

Finally, if `H0 y0 = T` is the landed positive reference solve and `H=H0+E`,
then

`y-y0 = -H^-1 E y0`,

and

`||E||_2 <= ||E||_F <= delta sqrt(s_-^2+s_+^2)`.

Therefore every component of the actual squared-amplitude vector is bounded
below by

`min(y0) - ||H^-1||_2 ||E||_2 ||y0||_2`.

`CovariancePerturbationCertificate` fails closed unless this lower bound is
strictly positive. The concrete 2x2 solve is then evaluated only as the actual
coefficient value and is independently cross-checked in tests against
`numpy.linalg.solve`; positivity is not inferred merely from that numerical
solve.

The scalar upper/lower arithmetic uses outward `nextafter` steps after each
positive product/sum/division or lower-bound subtraction. This is a conservative
binary64 implementation boundary, not an interval-arithmetic or Lean proof.

## Official Lean cross-check boundary

Repository: `openai/NavierStokesAndEuler`

Pinned commit:
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

`NavierStokes/Covariance.lean` formalizes the exact signed reference matrix,
determinant, inverse formula, strict-cone positivity and square-root amplitudes.
Its module header explicitly states that the actual integrated columns include
approximation errors and that the file does **not** identify those columns with
the exact model or prove the analytic error estimates. Accordingly, this Python
perturbation layer is a separate finite-dimensional verifier, not a claim that
the pinned Lean already proves Eq. (7.28) for the actual pulses.

## Tests and deliberate boundary

`tests/test_stress_cone_perturbation.py` uses synthetic perturbation fixtures
only. It independently checks the determinant with `numpy.linalg.det`, the
spectral inverse norm with `numpy.linalg.norm(inv(H), 2)`, and the coefficients
with `numpy.linalg.solve`. It also checks failure when the supplied concrete
column errors exceed the certified bound, when the determinant margin is lost,
and when determinant nondegeneracy survives but the conservative positive-cone
margin does not.

This increment still does **not**:

1. construct the actual Eq. (7.27) pulse integrals from the paper-exact
   Proposition 5.5 background, phases, damping and homogeneous pulses;
2. prove the uniform analytic `|e_sigma| <= C S_*^(-1/2)` estimate of Eq. (7.28);
3. derive the paper's full differentiated Eq. (7.29) bounds or shell-edge
   estimates;
4. construct the paper-exact target `T_{0,*}`;
5. solve the amplitude transport/ODE, realize the oscillation by a supported
   curl, or perform mean correction / Section 9 residual iteration.

Thus caller-supplied error vectors and bounds remain certified inputs to a
**formal-structure** implication. They are never promoted to paper-exact data.
