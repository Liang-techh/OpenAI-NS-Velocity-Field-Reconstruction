# Scientific status

The release contains two complete finite-window research candidates, not accepted Navier-Stokes solutions. ST054-Q2 is the lower-volume-L2 default among the two included ST054 candidates. ST054-M3 is the lower-peak alternative. Neither dominates every metric. Newer ST055 comments did not supply a complete candidate-and-evidence package at selection time and are not promoted.

## Fixed problem

The full residual is `R = u_t + (u dot grad)u + grad(p) - 0.01 Delta(u) - f`.
Viscosity is 0.01; time is [0.25,0.75]; the physical domain is R^3 and the evaluation box is [-2,2]^3. Velocity and pressure have smooth compact support inside r<2 and |z|<2. The independently prescribed two-parameter curl force has a,c in [0,10]. Initial energy E(0.25)=1 prevents trivial amplitude collapse.

Both momentum gates require sampled full-vector max and fixed-time spatial volume L2 below 0.001. L2 is `sqrt(64 * mean(sum(R**2, axis=1)))`, not component RMS, training loss or a space-time average. Divergence, energy, boundary, core-sign and one-probe scaled-core gates remain those of runtime/validate.py.

Both candidates fail the momentum gates. Sampled maxima are not continuous upper bounds; finite checks do not certify exact source identity or blow-up. Source-paper annular pulses and matched heat exterior are not established in this family. Small core direction tests do not establish whole-field geometry.

## Evidence boundaries

See evidence/source_results.json for paired seeds, domain and original scope. ST054 changed pressure only relative to each parent: its velocity geometry did not improve as a consequence of that fit. Core drift and support-edge peaks remain. Full independent replays in evidence/replay are executed during the release build, without refitting, using the original samples and operator.

Scientific failures return exit 1. A green CI job deliberately checks that the rejection remains reproducible. It does not rename failure as success. No paper-exact, OpenAI-field identity, singularity, global optimum or continuous residual certificate is claimed.
