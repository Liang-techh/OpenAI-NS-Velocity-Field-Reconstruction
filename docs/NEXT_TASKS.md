# Next constructive tasks and acceptance criteria

The next priority is **actual profile data**, not more tests of the toy exponent.
Preserve the existing numerical regression suite while implementing these dependencies.

## P0: Instantiate one admissible leading-profile construction

Translate the constructive choices behind Theorem 4.6 and Appendices A/B/C.
Produce the outer heat profile, inner analytic profile, matching moments,
shear modification and cone constraints. The pointwise moment solver is a
reusable subroutine, not a substitute for that construction.

Deliver a deterministic constructor and a manifest containing every free
parameter, equation reference, cutoff, support radius, quadrature order,
coefficient/table hash and truncation tolerance. Save the actual moment
systems B(eta), Q_eta, d(eta) and account for invertibility and smallness
uniformly on the required eta domain. Establish parity, axis regularity,
pressure balance, matching and support requirements. Do not set
`paper_exact=True` solely because a numerical plot or sampled constraints
look right.

## P1: Solve the background recursion from that profile

Implement Eqs. (5.2)-(5.6) and the compact moment corrections. Use the
existing analytic cutoff/curl assembly only after generating its genuine
coefficients. Record why the cutoff schedule is admissible. Compare
residuals and derivatives while increasing both coefficient order and
numerical resolution; separate truncation, quadrature and roundoff errors.

## P2: Instantiate Sections 6-9

Implement the dyadic charts, phase transport, stress cone and amplitudes,
curl remainders, compact mean corrections, and residual-improvement cycle.
Each stage needs independently checked conservation, support, moment and
residual identities. A manually tuned wave-frequency list or finite plot
cannot establish the required iteration and all-order convergence.

## P3: Final compact field, pressure and smooth force

Assemble Eq. (10.4) using the actual local field. Supply the spatial/time
cutoffs and independently establish the force extension through t=1.
Check the paper's support, energy and singular-path requirements. Only then
replace the top-level incomplete status with an evidence-backed completion
state. Numerical residuals supplement the analytic/formal evidence; they do
not replace it.

## Upstream verification and handoff

Pin and review the relevant official Lean source definitions and theorem
statements; record source paths, theorem names and hashes. Run the actual
Lean build before claiming a successful formal cross-check. The source pin
currently present is observed commit metadata only.

Existing edits must be integrated on a branch and tested against the newest
remote revision. This delivery is a patch against refreshed snapshot
`b01aaeebc6ec8a39f6692109b6a3c81b92cd0f32`; it does not overwrite or supersede
subsequent work by other agents. Never use forced reset/push to apply it.
