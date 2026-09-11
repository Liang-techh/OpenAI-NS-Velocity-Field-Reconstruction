# Section 7 signed stress-cone reference provenance

Status: **formal-structure**.  This file does not upgrade Stage 3/4 or the
repository-wide paper-exact gate.

## Published construction mapped here

Paper: OpenAI, *Finite Time Blowup for Navier-Stokes* (September 2026),
Section 7, especially Proposition 7.5 and Eqs. (7.24)-(7.30) (printed pp.
82-83 in the inspected PDF).

The paper first forms, for each slow box and signed pair, the two-column
pulse-covariance matrix `H`.  Its normalized columns are compared with a
reference signed pair, the order-zero target is placed strictly inside their
cone, and the two positive squared amplitudes are obtained by applying the
inverse matrix.  The amplitudes are then square-rooted in Eq. (7.24), and their
covariance realizes the target to the stated scale in Eq. (7.26).  The actual
paper proof also controls the difference between the integrated columns and the
reference columns and derives uniform inverse/amplitude derivative bounds.

## Pinned official Lean cross-check

Repository: `openai/NavierStokesAndEuler`

Pinned commit:
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

Primary module used by this increment:
`NavierStokes/Covariance.lean`.

The runtime module `src/openai_ns_reconstruction/stress_cone.py` transcribes the
finite-dimensional reference algebra from that file:

- signed columns `(-a,-b) scaleMinus` and `(-a,+b) scalePlus`;
- target `(-m,t)`;
- determinant `-2 a b scaleMinus scalePlus`;
- squared coefficients
  `(b m-a t)/(2 a b scaleMinus)` and
  `(b m+a t)/(2 a b scalePlus)`;
- strict positivity iff `|a t| < b m` under positive column parameters;
- amplitudes as the positive square roots of those coefficients;
- the scaled covariance identity produced by multiplying amplitudes by
  `sqrt(epsilon) * mask`;
- the manuscript-ratio implication using
  `a=-c_* sqrt(1+u_*^2)` and
  `|c_* t/m| < u_*/sqrt(1+u_*^2)`.

`tests/test_stress_cone.py` independently compares the explicit coefficient
formula with `numpy.linalg.solve`, compares the determinant with
`numpy.linalg.det`, reconstructs the scaled covariance by direct matrix
multiplication, and checks fail-closed behavior at and outside the strict cone.
The independent solve is a numerical cross-check only, not a theorem proof.

## Deliberate boundary

This increment is **not** the full Proposition 7.5 stress realization.  In
particular it does not:

1. identify the paper's actual pulse-integrated covariance matrix `H` with the
   ideal signed matrix;
2. construct `H` from the true Proposition 5.5 background, slow boxes, phases,
   damping and homogeneous pulses;
3. prove the Eq. (7.28) approximation of the two normalized columns or the
   determinant/inverse and positive-amplitude bounds of Eq. (7.29);
4. construct the paper-exact order-zero target `T_{0,*}` from the residual;
5. prove parameter/slow-variable derivative bounds or shell-edge vanishing;
6. solve the later pulse amplitude ODE, realize the physical oscillation by a
   supported curl, or perform compact mean correction / Section 9 iteration.

Therefore caller-supplied `a`, `b`, scales, target components, `c_*`, `u_*`,
`epsilon`, or masks remain formal inputs.  Passing the strict cone test does
not make them paper-exact data, and no surrogate wave is produced by this
module.
