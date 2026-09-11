# Section 7 pulse covariance error-budget provenance

Status: **formal-structure**. This increment does not upgrade Stage 3/4 or the
repository-wide paper-exact gate.

## Published construction mapped here

Paper: OpenAI, *Finite Time Blowup for Navier-Stokes* (September 2026),
Section 7, especially Eqs. (7.2), (7.21), (7.27), (7.28), and Proposition 7.5.

The signed phase-shear variable is affine in the pulse coordinate,

` s(v) = sigma (u_*/2 + u_* v/L_s) `,

so at the pulse center `v=L_s/2`, `s=sigma u_*`.  Lemma 7.4 writes the
homogeneous tangential pulse as

` t^h_sigma = x (e_r - s(v) K) + y N `

and gives

` y/x = c_* sqrt(1+s(v)^2) + O(S_*^-1) `.

After the positive factors in Eq. (7.27) are divided into the scalar `h_sigma`,
the normalized covariance column is therefore a positive weighted average of

` g(s(v)) + r(v), `

where, in frozen `(N,K)` coordinates,

` g(s) = (c_* sqrt(1+s^2), -s) `

and `r` is the ratio/frame error.  At the center,

` g(sigma u_*) = (-A_c, -sigma u_*), `

with `A_c=-c_*sqrt(1+u_*^2)>0`, exactly the reference vector in Eq. (7.28).

## Executable implication

`SignedPulseCovarianceErrorBudget` accepts two quantities that must already be
certified analytically upstream:

1. `ratio_error_bound`, satisfying `|r(v)| <= E_ratio` on the full pulse
   support;
2. `normalized_first_moment_bound`, satisfying
   `E_w[|v-L_s/2|/L_s] <= M_1` for the positive Eq. (7.27) weight
   `w(v)=psi(v)^2 x(v)^2`.

The derivative of the ideal direction obeys

` |g'(s)|^2 = 1 + c_*^2 s^2/(1+s^2) <= 1+c_*^2. `

Because

` |s(v)-sigma u_*| = u_* |v-L_s/2|/L_s, `

the mean-value theorem and positivity of the normalized weight give

` |E_w[g(s(v))]-g(sigma u_*)|
    <= u_* sqrt(1+c_*^2) M_1. `

Adding the averaged ratio defect yields the executable bound

` |e_sigma|
    <= E_ratio + u_* sqrt(1+c_*^2) M_1. `

The implementation outward-rounds the positive scalar arithmetic in binary64
and preserves exact zero contributions.  `error_vector_for` is only a
consistency check for a supplied normalized column; it does not establish that
the supplied vector is the continuous Eq. (7.27) pulse integral.

For both signs, `uniform_two_sign_error_bound` takes the maximum of the two
analytic budgets after checking that they use the same `c_*` and `u_*`.  That
single number is the `delta` expected by the already-landed
`CovariancePerturbationCertificate`; the concrete error-vector ordering in that
module follows transverse sign, whereas the manuscript label `sigma` satisfies
`(-A_c,-sigma u_*)`, so callers must not identify those names by string alone.

## Recovering the Eq. (7.28) square-root rate

The paper obtains an `O(S_*^-1)` homogeneous-ratio error from Lemma 7.4 and an
`O(S_*^-1/2)` normalized first moment from the Gaussian pulse concentration in
the proof of Proposition 7.5.  If

` E_ratio <= C_ratio/S_* `

and

` M_1 <= C_moment/sqrt(S_*) `,

then for `S_*>=1`,

` |e_sigma|
   <= [C_ratio + u_*sqrt(1+c_*^2) C_moment]/sqrt(S_*). `

`sqrt_scale_error_envelope` implements precisely this scalar rate collapse.  It
does not prove either component estimate for the actual pulse.

## Tests and independent cross-check boundary

`tests/test_pulse_covariance_budget.py` constructs only synthetic positive
weights and tangential-ratio values.  It independently forms the weighted
average with NumPy, splits the actual synthetic error into averaged local defect
plus center drift, and checks that each piece lies under the analytic budget.
It separately samples the closed-form derivative norm of `g` to cross-check the
global Lipschitz constant and verifies the `C/sqrt(S_*)` scalar envelope.

These fixtures test the implication and sign conventions only.  They are not
surrogate waves and are not paper pulse data.

## Remaining paper-exact boundary

This increment still does **not**:

1. instantiate the true Proposition 5.5 background on each slow box;
2. construct the actual homogeneous pulse functions `x,y,psi` or prove their
   Lemma 7.4 uniform ratio/frame estimate;
3. prove the Gaussian weighted first-moment estimate for those actual pulses;
4. evaluate the continuous Eq. (7.27) columns `H_sigma` and their positive
   scales `h_sigma` from paper-exact data;
5. identify caller-supplied vectors with those columns;
6. construct the paper-exact target `T_{0,*}`, amplitude ODE, supported-curl
   oscillations, mean corrections, or the Section 9 residual iteration.

Thus the new module closes only the theorem-shaped analytic implication that
turns certified Lemma-7.4 and pulse-concentration inputs into the Eq. (7.28)
normalized error budget.  Stage 3/4 remains **formal-structure** and
`paper_exact_velocity_available` remains false.
