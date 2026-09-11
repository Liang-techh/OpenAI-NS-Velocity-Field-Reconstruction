# Section 7 cylindrical coefficient-curl jet provenance

Status: **formal-structure**. This increment does not upgrade Stage 3–6 or the repository-wide paper-exact gate.

## Construction mapped here

Paper: OpenAI, *Finite Time Blowup for Navier–Stokes* (September 2026), Section 7, at the vector-potential / curl realization step following the primary oscillatory amplitude construction.

Pinned official formalization repository: `openai/NavierStokesAndEuler`

Pinned commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Primary modules cross-checked:

- `NavierStokes/CurlClassBounds.lean`, especially Formula (30), `normalCoefficient`, `inverseCarrier`, `curlRemainder`, `vectorPotential`, `realizedCoefficient`, and `cylindricalCurl_vectorPotential`;
- `NavierStokes/LocalizedCurlRealization.lean`, especially `RawData.native_realizes_curl`, `RawData.native_divergence_zero`, and the local/common zero-germ glue conditions.

The already-landed `curl_realization_algebra.py` records the finite-dimensional identity

`B = |n|^-2 (n × a)`,

with tangent realization

`a_corrected = a + (i/K) curl(B)`.

That layer intentionally accepted `curl(B)` as an arbitrary caller-supplied complex three-vector. The present module `src/openai_ns_reconstruction/cylindrical_curl_jet.py` narrows that interface. Callers now provide the coordinate derivative jet of the cylindrical coefficient components `(B_r,B_theta,B_z)`, and the module constructs

- `curl(B)_r = (1/r) ∂_theta B_z - ∂_z B_theta`,
- `curl(B)_theta = ∂_z B_r - ∂_r B_z`,
- `curl(B)_z = ∂_r B_theta + B_theta/r - (1/r) ∂_theta B_r`.

The connection term `B_theta/r` is kept explicitly. This is the genuine standard cylindrical curl of component data at `r>0`; treating the components as Cartesian would omit that term.

`TangentCylindricalCurlJet` also removes caller choice of the coefficient value itself: it reuses the exact-tangency gate from `TangentCurlRealization` and takes the value of `B` from the pinned `|n|^-2(n×a)` construction. Only the three derivative vectors remain upstream inputs.

## Independent regression check

`tests/test_cylindrical_curl_jet.py` validates the cylindrical formula against an independent coordinate route. It defines an explicit smooth complex cylindrical vector field, converts it to Cartesian components, differentiates the Cartesian field by centered finite differences in `(x,y,z)`, forms the Cartesian curl, and compares that result with the production cylindrical formula transformed back to Cartesian coordinates.

This check is intentionally not production-formula-versus-production-formula. Separate tests verify the tangent Formula (30) assembly and fail-closed behavior for `r<=0` and non-tangent amplitudes.

## Deliberate boundary

This increment still does **not** construct a paper-exact divergence-free oscillatory wave. In particular it does not prove that the supplied derivative jet is the derivative jet of the actual localized paper coefficient. It also does not establish:

1. the Proposition 5.5 paper-exact background/base field;
2. the actual Section 7 phase/frame hypotheses and pulse-amplitude datum on every active slow box;
3. analytic differentiation of the genuine localized product `cutoff * amplitude` into the derivative jet accepted here;
4. smoothness and nonvanishing of the true phase normal on the native phase patch;
5. support / zero-germ hypotheses needed by `LocalizedCurlRealization.RawData`;
6. the common-potential glue across cells;
7. the theorem-level `div(curl)=0` conclusion for the assembled physical wave;
8. compact mean corrections or the Section 9 residual-improvement iteration.

Therefore this is reusable derivative infrastructure only. `paper_exact_velocity_available` must remain `false` until the upstream derivative data and the remaining support/gluing hypotheses are instantiated from the actual paper construction.
