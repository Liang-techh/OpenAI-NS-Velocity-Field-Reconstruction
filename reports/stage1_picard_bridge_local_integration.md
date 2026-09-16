# Stage 1 finite Picard bridge local integration report

Status: **provisional formal-structure increment; focused acceptance passed**.

This report records the finite Picard family bridge from the actual solver
anchor. It does not claim a converged profile, global weighted coefficient
state, final pressure, or a fixed point.

## Repository identity

- Baseline main commit: `67ab5ecd32f80225df8bce2760b943c20a7d3d7b`
- Local integration branch: `codex/stage1-complete-slow2`

## API and scope

`formal_axis_picard_family_state(solver)` returns
`FormalAxisPicardFamilyState`. Its
`jet_prefix(iterations, max_n, m, eta)` returns finite angular/axial sparse
mixed-scale rows `0..max_n`; `iterations=0` is the actual reference `x0`.
Each later step applies `reference + R(input)/(2 Lambda)` through the same
pinned graph, with full eta-derivative family semantics. `jet_pair` selects
one row. The solver-level `picard_map_families` and
`picard_map_jet_prefix` expose the one-step routed map for independent use.

The bridge keeps row and eta caches local to each call. It has no fixed
iteration cutoff and no persistent unbounded family cache. Updated-iterate
pressure is intentionally outside this API; the private graph pressure value
under an external input pair is `P(input)`.

The wide physical assembly now exposes the sparse maps `F`, `E`, `U`,
`dU_deta`, `average_U`, `d_average_U_deta`, `V0`, and `Pi`. Its conditional
tail mapping uses exact Fraction scaling for inverse-Lambda fields and an
upward-verified `sqrt(2 X)` factor for `E`; it excludes coefficient and
axis-pressure quadrature roundoff.

## Evidence and proof boundary

The focused bridge acceptance command was:

```text
python -m pytest -q tests\test_axis_coefficient_picard_family.py tests\test_axis_coefficient_formal_solver.py tests\test_axis_coefficient_profile_prefix.py tests\test_natural_axis_wide.py -W error
```

It exited `0` with `13 passed` in `25.79 s`. The target independently covers
finite iterations through row `2`, x1/x2 agreement, the formal profile prefix,
wide `d_average_U_deta`/`V0`, and exact conditional tail scaling. No full-suite
rerun was performed.

The source-pinned manual proof in
`references/AXIS_COEFFICIENT_PICARD_FAMILY_PROVENANCE.md` is accepted. Under
the theorem-side compatible fixed-point hypotheses, it establishes that formal
row `n` stabilizes after update `n` and equals the compatible fixed-point
coefficient. This is an exact-arithmetic filtration result only; it does not
bound errors in the evaluated scalar jets, Decimal rounding, or coefficient
roundoff. The jet values are finite, but their numerical errors lack certified
bounds. Evaluated scalar jet enclosures and a coefficient-roundoff budget remain
the next Stage 1 blocker.

Global compatible `AxisSpace` membership, adjacent-eta error control,
convergence, a complete physical profile, and paper-exact status remain open.
