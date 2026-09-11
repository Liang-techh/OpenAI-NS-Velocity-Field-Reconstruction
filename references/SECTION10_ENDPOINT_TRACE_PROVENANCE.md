# Section 10 endpoint-trace enclosure provenance

Status: **formal-structure only**. This increment does not make the reconstructed velocity or force paper-exact.

## Paper / Lean boundary

This adapter sits between the already-landed late-time derivative-majorant infrastructure and the endpoint jets required by the Section 10 Taylor--Borel glue. The relevant pinned formal interfaces are `CandidateFromLimits.lean`, which assumes locally uniform limits of every full spacetime derivative of the actual closed-past Navier--Stokes residual, and `SpacetimeGluing.lean` / `SpatialBorelExtension.lean`, which consume the resulting endpoint jets.

The repository already records a sufficient conditional hypothesis

`sup_{x in K_m} ||partial_t D^n R(t,x)|| <= C (1-t)^(-alpha)`, with `alpha < 1`,

and integrates it on the official `t >= 3/4`, `T = 1` plateau. `section10_endpoint_trace.py` takes no new fitted constants. For each finite derivative degree in a certified `Section10EndpointMajorantLadder`, it evaluates the supplied pre-endpoint full residual jet at one time and returns the closed dense-tensor ball with radius equal to the theorem-side endpoint-tail budget. Conditional on the independently proved majorant, the true endpoint limit must lie in that ball.

## Executable convention

Full order-`n` spacetime derivatives use the existing dense shape `(3,) + (4,)*n`, with derivative coordinates ordered `(time,x,y,z)`. Enclosure distances use the Euclidean norm of this dense coordinate tensor. Any upstream majorant used with this adapter therefore has to be certified in that norm or in a norm that dominates it.

The regression oracle is independent of the implementation: it uses an explicit analytic family `L_n + a_n (1-t)^beta_n e_n`, whose endpoint value and exact derivative norm are available in closed form. The test supplies a strict 1% analytic upper bound on that derivative norm, so the expected enclosure radius is `1.01 |a_n| (1-t)^beta_n`; the exact endpoint remainder is checked independently against `|a_n| (1-t)^beta_n`.

## What this closes

- prevents the finite-order endpoint-majorant ladder from remaining only scalar bookkeeping;
- provides a constructive, shrinking approximation/enclosure for each full endpoint residual jet through the ladder's finite degree;
- reuses the official `t>=3/4`, `T=1` gate and the existing full-spacetime tensor validator;
- performs no sampling-based fitting or majorant shrinking.

## What remains open

- derive the majorants from the actual completed Section 9/10 residual rather than caller-supplied hypotheses;
- prove the queried spatial point belongs to the compact window represented by `spatial_window` (this adapter carries the index but does not reconstruct that set membership);
- pass from finite-degree enclosures to actual locally uniform endpoint-limit functions at every order;
- derive genuine analytic `templateBound` values for the Borel scale schedule from those endpoint jets;
- prove the endpoint trace/normal-jet matching, all-order smooth future glue, compactly supported smooth forcing through `t=1`, blow-up path and uniform bounded energy.

Accordingly `paper_exact_velocity_available` must remain `false`.
