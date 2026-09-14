# Section 9 finite-prefix field assembly provenance

Status: **formal-structure only**.

This increment makes the already-admitted finite-stage ledger executable once
actual stage field callables are supplied.  It does not construct those missing
paper stages and therefore does not promote the reconstruction to paper-exact.

## Pinned formal source

Repository: `openai/NavierStokesAndEuler`

Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Files and definitions:

- `NavierStokes/DiagonalJetBounds.lean`
  - `DiagonalJetBounds.uncutPrefix A N x = sum j in range N, A j x`.
- `NavierStokes/MixedDiagonalResidual.lean`
  - `MixedDiagonalResidual.uncutVelocity A B J = curl(uncutPrefix A (J+1)) + uncutPrefix B (J+1)`.
- `NavierStokes/WholeDomainActualStageBounds.lean`
  - `WholeDomainStageBounds.finiteVelocity B N0 hN J` is the above uncut velocity for the actual `potentialStages/directStages`.
  - `WholeDomainStageBounds.finitePressure B N0 hN J = uncutPrefix pressureStages (J+1)`.

The Python constructor therefore requires the materialized stages to be exactly
`0,...,J`, never silently inserting zero for a missing stage.  For each stage it
consumes the spatial curl of the potential, the direct velocity and the pressure
as callable fields.  The spatial curl is deliberately supplied by the stage
producer rather than reconstructed numerically from samples.  Linearity then
gives the literal finite-prefix velocity algebra.

The stage family identities must agree with an existing
`Section9ActualStageEstimatesAdmission`, tying the executable prefix to the
canonical finite velocity/pressure ledgers already admitted from
`ActualStageEstimates.stageEstimates_of_representations`.

## Deliberate non-claims

The repository still does not have the actual paper stage callables produced by
the full Sections 6--8 correction cycle.  The test callables are explicitly
fixtures and are not evidence for those fields.  This increment therefore does
not claim:

- that the actual stage fields are materialized;
- that the local rank increment or post-rank remainder values are materialized;
- that the finite prefix has converged as `J -> infinity`;
- that Eq. (9.21)'s infinite/local limiting field is closed;
- that the residual is small merely because a finite prefix can be evaluated;
- or that paper-exact velocity is available.

The next blocker is to make Agent 3/Section 8 output genuine materialized stage
fields with the canonical ledger identities, then feed those fields through this
prefix constructor and the existing stage-estimate/selected-schedule chain.  A
separate convergence argument remains mandatory before Agent 4 may consume an
actual limiting local field.
