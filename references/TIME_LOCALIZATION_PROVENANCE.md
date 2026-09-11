# Section 10 time-localization provenance

Pinned formal source: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Relevant modules:

- `NavierStokes/SmoothCutoffs.lean`: defines `timeSwitch t = 1 - scaledCutoff (4/3) t`; proves `timeSwitch=0` on `|t|<=3/8`, `timeSwitch=1` for `t>=3/4`, and `0<=timeSwitch<=1`.
- `NavierStokes/TimeLocalization.lean`: defines activated velocity/pressure by multiplying both by `timeSwitch`; proves early zero, late identity, preservation of spatial divergence-free structure, preservation of terminal speed blowup, and the exact residual identity
  `R(chi u, chi p) = chi R(u,p) + chi' u + (chi^2-chi) (u.grad)u` on the presingular interval.

Executable artifact: `src/openai_ns_reconstruction/time_localization.py`.

Truth boundary:

- The executable switch uses this repository's explicit C-infinity even bump representative with the same inner/outer radii as the pinned Mathlib `ContDiffBump`. Hence the zero region `|t|<=3/8` and the unit region `|t|>=3/4` are exact for the executable representative.
- We do **not** claim pointwise equality in the transition collar with Mathlib's noncomputable `ContDiffBump`.
- `time_switch_derivative` is the exact analytic derivative of the executable representative.
- `activated_residual_formula_numeric` evaluates the pinned theorem's right-hand side using analytic `chi,chi'` and independent numerical derivatives of the unactivated fields. Tests compare it against a separate numerical residual evaluation of the activated fields; the expected value is not defined as `R-f` after setting `f=R`.
- This increment does **not** construct the unresolved incoming local fields and does **not** construct the separate smooth global force extension through `t=1`. It therefore does not make the final velocity paper-exact.
