# Physical-chart finite prefactor provenance

This increment belongs to the Issue #2 / Agent 7 all-order analytic-closure lane and is stacked directly on PR #345.

Pinned source: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, `NavierStokes/SlowBorelBase.lean`.

The formal chain used here is exact:

- `physicalChart_jet_bound` supplies, for every fixed derivative order `i`, an existential positive constant `C_i` on a fixed bounded `X` interval.
- `physicalChart_finite_bound` defines `D = 1 + sum_{i=0}^m C_i` and proves all physical-chart jets through order `m` are bounded by `(D/q)^i`.
- `physical_composition_bound` contributes the exact chain-rule factor `m!`.
- `exists_physical_uncut_tail` combines that with the ordinary tail factor `2^-J q^(P+m)` to obtain the physical bound with finite prefactor `m! * 2^-J * D^m`.

The executable certificate deliberately does **not** invent numerical values for the existential `C_i` or `D`. It reconstructs their symbolic dependency graph from the pinned theorem, derives the `X` interval from the same hierarchy-owned common-support certificate used by the `C[j,m]`/recurrence/cutoff chain, and binds the exact `m!`, `2^-J`, `D^m`, first-omitted order, and physical q-power into one finite request.

No caller can supply `C_i`, `D`, a replacement support interval, a second cutoff schedule, a first-omitted exponent, or an independent `C[j,m]` table through this API.

Truth boundary:

- finite physical-chart bound **shape** replayed from the pinned theorem: yes;
- physical composition prefactor **shape** replayed exactly: yes;
- numerical existential physical-chart constant materialized: no;
- Lean build/proof replay in Python: no;
- actual PDE residual evaluated: no;
- total all-order hierarchy: no;
- infinite diagonal schedule: no;
- all-jets flatness / super-algebraic residual convergence: no;
- paper exactness / full reconstruction: no.

Upstream return to Agent 2 is unchanged and sharper after PR #348: closing the finite Picard derivative triangle is not a total coefficient solver. The eventual backend must expose, for every requested positive order, one actual assembled/repaired coefficient object together with theorem-backed derivative majorants, common compact support, and exact retained recurrence evidence under one coefficient state. Once that provider is total, this symbolic physical-chart constant is sufficient for the theorem's finite-per-request prefactor; no uniform numerical `D` over all derivative orders is required by the pinned quantifiers.
