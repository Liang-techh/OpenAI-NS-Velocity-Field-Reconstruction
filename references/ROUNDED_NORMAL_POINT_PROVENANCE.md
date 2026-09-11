# Rounded-normal point adapter provenance

Source pin: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, `NavierStokes/PhaseEstimates.lean`.

This increment transcribes the constructive input/output boundary around the pinned theorem `rounded_normal_estimates` and the downstream smallness gate used by `actual_phase_geometry`.

Implemented in `src/openai_ns_reconstruction/phase_box_certificate.py`:

- `lean_nonzero_round` / `lean_rounded_frequency` reproduce the official Lean punctured-lattice rounding used by the quantitative theorem: floor `k*target`, except a zero floor is replaced by `1`.
- `RepresentativePhaseData` derives `g=(R0*FR0,GR0)`, an orthogonal unit `K`, and the representative frequency `B*(K-(sigma*u/(L*||g||^2))*g)`. Therefore the theorem's structural equalities `g=[R0*FR0,GR0]`, `K ⟂ g`, and `[target/R0,pz]=representativeFrequency ...` are satisfied by construction rather than accepted under a floating equality tolerance.
- `RoundedNormalPointCertificate` fails closed unless the numerical inequalities in `rounded_normal_estimates` hold at the supplied point: scale/domain bounds, frequency/component bounds, inverse-radius/shear bounds, radius diameter, and the `FR/GR` representative derivative error budget `M*(S^-3+epsilon^2)`.
- The adapter independently evaluates `explicitNormal`, `normalVelocity`, the theorem envelope `phaseConstant(M)*phaseError(S,epsilon,k)`, and the lower/frame quantities used after the additional `delta<=B/2` gate.

The module is intentionally **not** a uniform slow-box certificate. It does not infer local C1/C2 bounds from samples, nor does it claim that caller-supplied derivatives arise from the paper-exact base field. Uniform Eqs. (7.9)-(7.11) still require an analytic/local-base layer proving the `LocalBaseBounds` hypotheses over every active box and connecting those functions to the unresolved paper-exact leading/background construction.

This artifact remains `formal-structure`; it does not create a paper-exact oscillatory wave, stress decomposition, amplitude solve, mean correction, or Section 9 residual iteration.
