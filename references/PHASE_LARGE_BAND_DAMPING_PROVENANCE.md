# Section 7 large-band damping-error admission provenance

Status: **formal-structure**

This layer records a fail-closed application identity for the pinned Lean theorem
`BasePhaseGeometry.FamilyData.damping_error` after the existing large-band
signed-family/base-source/scalar admission has succeeded.

## Pinned formal source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Module: `NavierStokes/BasePhaseGeometry.lean`
- Theorem: `BasePhaseGeometry.FamilyData.damping_error`

At the pinned source, the theorem concludes on an active carrier that

`|nu * rbar^-2 * viscosityDamping(i)
 - referenceViscosity(M,r0Center)^-1 * viscosityDamping(i)|
 <= dampingConstant(M) / D.scale(i)`.

The same source defines

`dampingConstant(M) = 4*M*(6*M+5)*phaseConstant(M)`.

The repository's existing `phase_frame_bounds.py` / `phase_large_band_frame.py`
already transcribe the scalar specialization and map `D.scale(i)` to the
large-band scale `S_*=ell^2`.  The new admission layer does not alter those
constants; it attaches theorem-level identity/provenance requirements to the
actual damping expressions.

## Accepted evidence

Only `analytic-theorem` and `formal-theorem` witnesses are accepted.  Sampled,
fitted, or numeric-scan evidence is rejected.

A successful witness must certify all of the following for the same admitted
signed family/source revision:

- active-carrier membership used by the theorem;
- identity of the actual term `nu*rbar^-2*viscosityDamping(i)`;
- identity of the reference-viscosity damping term;
- identity of the common `viscosityDamping(i)` family quantity;
- reference-viscosity and `r0Center` identities;
- `D.scale(i)=S_*=ell^2`;
- dependency on the already-admitted family input;
- exact pinned repository, commit, and theorem name.

The implementation additionally checks that the existing scalar damping
envelope agrees with the independently named
`dampingConstant(M)/S_*` specialization.

## Explicit non-claims

This layer does **not** materialize an actual damping field, prove the genuine
Section 6 carrier/slot family, provide the Proposition 5.5 paper-exact
background, or verify the theorem for every active family index.  It therefore
does not establish Eqs. (7.9)-(7.11), construct a stress-cone amplitude/wave,
or make a paper-exact velocity available.

Truth flags remain false for:

- `actual_damping_field_materialized`
- `damping_error_machine_verified`
- `all_active_family_damping_verified`
- `uniform_eq_7_9_to_7_11_verified`
- `paper_exact_velocity_available`
