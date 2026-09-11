# Section 7 oscillatory curl-realization algebra provenance

Status: **formal-structure**.  This increment does not upgrade Stage 3/4 or the
repository-wide paper-exact gate.

## Published / formalized construction mapped here

Paper: OpenAI, *Finite Time Blowup for Navier-Stokes* (September 2026),
Section 7, at the vector-potential / curl realization step following the
primary oscillatory amplitude construction.

Pinned official formalization repository:
`openai/NavierStokesAndEuler`

Pinned commit:
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Primary modules cross-checked:

- `NavierStokes/CurlClassBounds.lean`, especially `normalCross`,
  `normalCoefficient`, `inverseCarrier`, `curlRemainder`, `vectorPotential`,
  `realizedCoefficient`, `cylindricalCurl_vectorPotential`, and
  `realizedCoefficient_divergence`;
- `NavierStokes/LocalizedCurlRealization.lean`, especially the preservation of
  tangency under the localized cutoff and the local/common curl/divergence
  identities once the genuine smooth cylindrical geometry is available.

The official Lean file calls the vector-potential identity **Formula (30)**.
For a real phase normal `n` and complex raw coefficient `a`, the executable
module `src/openai_ns_reconstruction/curl_realization_algebra.py` transcribes
only its coefficient-level algebra:

- `normalCross(n,a) = n x a`;
- `normalCoefficient(n,a) = |n|^-2 (n x a)`;
- the vector triple-product identity
  `n x normalCoefficient = n (n.a)/|n|^2 - a`;
- therefore the principal harmonic curl coefficient is the tangential
  projection `a - n(n.a)/|n|^2`, and it equals `a` only under exact tangency
  `n.a=0`;
- `inverseCarrier(K)=i/K` and `phaseFactor(K)=iK`, hence their product is `-1`
  for `K != 0`;
- the displayed coefficient-derivative remainder is `(i/K) curl(B)` and the
  tangent specialization of the realized coefficient is
  `a + (i/K) curl(B)`.

The production adapter deliberately requires exact `n.a == 0` before exposing
the tangent specialization.  It does not accept a numerical tolerance as a
substitute for the theorem hypothesis.  Approximate/non-tangent data may use
the projection routine, where the normal defect remains explicit.

## Independent regression checks

`tests/test_curl_realization_algebra.py` cross-checks the implementation against
independent NumPy dot/cross/projection arithmetic.  It also checks a complex
coefficient whose real and imaginary parts are both exactly tangent, verifies
`(i/K)(iK)=-1`, reconstructs the displayed remainder independently, and checks
fail-closed behavior for a zero phase normal, zero/nonfinite frequency,
nonfinite vectors, and non-tangent data.

These are finite-dimensional arithmetic regression checks.  They are not a
proof that any caller-supplied derivative is the paper's cylindrical curl.

## Deliberate boundary

This increment does **not** construct a paper-exact divergence-free wave.  In
particular it does not:

1. instantiate the phase normal and raw amplitude from the actual Proposition
   5.5 background / Proposition 7.5 pulse data;
2. solve the paper's pulse-amplitude ODE or construct the true order-zero target
   `T_{0,*}`;
3. construct or differentiate the actual localized coefficient `B` in the
   paper's cylindrical geometry;
4. establish smoothness, nonvanishing normal, support, zero-germ, phase-patch,
   or cutoff hypotheses required by `LocalizedCurlRealization.RawData`;
5. identify a caller-supplied `coefficient_curl` with the genuine cylindrical
   curl of `B`;
6. promote the algebraic split to the formal divergence-free conclusion, which
   in Lean follows only after the genuine curl identity and `div(curl)=0`
   hypotheses are satisfied;
7. perform compact mean correction or the Section 9 residual-improvement
   iteration.

Thus the module is reusable infrastructure for the future genuine wave
construction, but `paper_exact_velocity_available` must remain `false`.
