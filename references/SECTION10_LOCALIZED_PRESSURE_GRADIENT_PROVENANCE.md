# Section 9 -> Section 10 localized pressure-gradient provenance

Status: **formal-structure / bounded derivative adapter only**.

## Pinned paper/formal source

At the repository-pinned OpenAI formalization commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`,
`NavierStokes/SpatialLocalization.lean` defines

```text
cutPressure (p) := spatialCutoff * p
```

using the same fixed spatial cutoff that localizes the velocity potential.
The theorem file also proves support and smoothness properties for this cut
pressure and later feeds `cutPressure p` into the localized Navier--Stokes
residual.

## This increment

`section9_section10_localized_pressure_gradient.py` connects that pinned
definition to an already-admitted finite Eq. (9.21) prefix jet. At one
caller-supplied spacetime point it evaluates

```text
p_cut = c p
grad(p_cut) = c grad(p) + p grad(c)
```

where `c` and `grad(c)` are taken atomically from the existing fixed Section 10
spatial-localization module. The caller cannot inject a generic cutoff or a
mismatched cutoff derivative.

Production code contains no finite differences. The regression independently
finite-differences the full scalar product `c(x) p(x)` in a point where both the
radial and axial cutoff factors are in their transition collars. Separate
checks verify that the plateau preserves `grad(p)` exactly, and that both the
cut pressure and its gradient are exactly zero outside the official fixed
support.

## Truth boundary

The input pressure jet is still the provider-supplied finite
`Section9FinitePrefixJetCertificate`. Its physical base-point binding is caller
metadata, and the transition-collar cutoff values are still the repository's
geometry-faithful executable C-infinity representative rather than a claim of
pointwise identity with Mathlib's noncomputable `ContDiffBump`.

This increment therefore does **not** establish any of the following:

- the actual Section 7/8 correction-field exporter;
- the infinite/locally-finite Eq. (9.21) field;
- all-order one-sided convergence or a smooth extension through `t=1`;
- a genuine Navier--Stokes residual artifact or forcing;
- endpoint residual closure;
- finite-energy or blow-up closure; or
- paper-exact reconstructed velocity.

`residual_artifact_ready=false` and `paper_exact_velocity_available=false`
remain mandatory.
