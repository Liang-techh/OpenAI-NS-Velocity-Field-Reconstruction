# Section 10 support / fixed-time energy certificate provenance

Status: **formal-structure / implication certificate**, not paper-exact velocity data.

Pinned formalization source: `openai/NavierStokesAndEuler` at commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, module
`NavierStokes/SpatialLocalization.lean`.

## Formal source mapped here

The pinned Lean module defines the closed support cylinder

- `supportCylinder := {x | radialSquare x <= 1/16 and |x 2| <= 1/4}`,

and proves, among other statements:

- `isCompact_supportCylinder`,
- `spatialCutoff_tsupport_subset`,
- `cutVelocity_tsupport`, and
- `cutVelocity_hasCompactSupport`.

Thus the spatially localized velocity has topological support inside a cylinder
of radius `1/4` and height `1/2`.  The Python geometry in
`spatial_localization.py` already mirrors those support constants.  The new
`section10_support_energy.py` adds an exact-rational guard for the same
geometry and records its elementary measure consequence:

`volume(cylinder) = pi * (1/4)^2 * (1/2) = pi/32`.

If, at one fixed time, an **independently certified** bound `|u| <= M` is
available on the actual localized field, then

`E(t) = (1/2) integral |u|^2 dx <= (1/2) * (pi/32) * M^2 = pi M^2 / 64`.

This implication uses only the support geometry and the supplied analytic
bound.  `tests/test_section10_support_energy.py` independently checks the
rational geometry identities and uses a labelled bounded smooth-cutoff fixture
with positive quadrature weights to cross-check the inequality numerically.
The fixture is not a surrogate for the paper field.

## Deliberate fail-closed boundary

This module does **not** infer `M` from grid samples and does not claim that an
arbitrary finite sampled maximum is a certificate.  In particular, it does not
prove the paper's stronger bounded-energy statement along the blow-up limit:
a pointwise `L-infinity` bound can itself diverge as `t -> 1-`, so the coarse
fixed-cylinder estimate may also diverge.  A genuine uniform energy result
still requires the actual completed local field and its paper-derived
space-time/scaling estimates.

It also does not certify endpoint smoothness, forcing closure, the blow-up
path, or the unresolved earlier profile/background/oscillatory construction.
Consequently `paper_exact_velocity_available` and `full_reconstruction` must
remain `false` until those upstream obligations are closed.
