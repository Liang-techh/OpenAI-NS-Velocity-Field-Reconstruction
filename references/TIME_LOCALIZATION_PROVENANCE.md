# Section 10 time-localization provenance

Pinned formal source: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Relevant modules:

- `NavierStokes/SmoothCutoffs.lean`: defines `timeSwitch t = 1 - scaledCutoff (4/3) t`; proves `timeSwitch=0` on `|t|<=3/8`, `timeSwitch=1` for `t>=3/4`, and `0<=timeSwitch<=1`.
- `NavierStokes/TimeLocalization.lean`: defines activated velocity/pressure by multiplying both by `timeSwitch`; proves early zero, late identity, preservation of spatial divergence-free structure, preservation of terminal speed blowup, and the exact residual identity
  `R(chi u, chi p) = chi R(u,p) + chi' u + (chi^2-chi) (u.grad)u` on the presingular interval.
- `NavierStokes/CandidateFromLimits.lean` / `NavierStokes/SpacetimeGluing.lean`: the later endpoint-force construction requires genuine locally uniform residual derivative limits as `t -> 1-`; the time-localization step must therefore be connected to that exact endpoint regime rather than to an arbitrary transition window.

Executable artifact: `src/openai_ns_reconstruction/time_localization.py`.

Truth boundary:

- The executable switch uses this repository's explicit C-infinity even bump representative with the same inner/outer radii as the pinned Mathlib `ContDiffBump`. Hence the zero region `|t|<=3/8` and the unit region `|t|>=3/4` are exact for the executable representative.
- We do **not** claim pointwise equality in the transition collar with Mathlib's noncomputable `ContDiffBump`.
- `time_switch_derivative` is the exact analytic derivative of the executable representative.
- `activated_residual_formula_numeric` evaluates the pinned theorem's right-hand side using analytic `chi,chi'` and independent numerical derivatives of the unactivated fields. Tests compare it against a separate numerical residual evaluation of the activated fields; the expected value is not defined as `R-f` after setting `f=R`.
- `section10_endpoint_localization_transfer(...)` is a fail-closed structural bridge from an already supplied `EndpointPowerLawMajorant` to the actual Section 10 endpoint window. It accepts only the pinned endpoint `T=1` and only validity windows starting at or after the exact unit plateau threshold `t=3/4`; transition-collar windows are rejected. On the certified interval it checks the executable correction coefficients exactly as `(chi, chi', chi^2-chi)=(1,0,0)`. Since activated and original fields are identical as functions throughout the late plateau, all residual spacetime derivatives transfer unchanged wherever they exist. This does **not** prove the source majorant, derive endpoint jets, or establish locally uniform convergence.
- The endpoint-transfer regression separately evaluates the Navier--Stokes residual on the original and activated fields at a late-plateau time, rather than defining a forcing by the same residual and checking `R-f=0`.
- This increment does **not** construct the unresolved incoming local fields and does **not** construct the separate smooth global force extension through `t=1`. It therefore does not make the final velocity paper-exact.
