# Pinned Section 5 uniform tail-order arithmetic provenance

Status: **formal-structure only**. Issue #2 / Agent 7 downstream all-order lane.

## Pinned source

Official formal repository: `openai/NavierStokesAndEuler` at
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Relevant definitions/theorems:

- `NavierStokes/SlowBorelBase.lean`
  - `exists_template_jet_bound`
  - `exists_admissibleScales`
- `NavierStokes/DiagonalScale.lean`
  - `exists_diagonal_scales`
  - `weighted_tail_bound`
  - `exists_uniform_tail_order`

`SlowBorelBase.exists_admissibleScales` obtains each normalized-template
constant from smoothness/compactness of the **actual coefficient at that stage**
and specializes the diagonal gain to `g(j)=2*h*j`. `DiagonalScale` then absorbs
each stage's finite family of `C[j,m]` values into its chosen cutoff scale and
takes the recursive doubling envelope.

This is an important quantifier clarification: the pinned scalar diagonal theorem
does not require one global numerical upper bound on all `C[j,m]`. It requires
the coefficient family to exist and be smooth at every order so that the
stagewise finite constants exist. The reconstruction still lacks that all-order
hierarchy provider; finite fixtures or a finite derivative frontier cannot fill
this quantifier.

## Increment

`background_uniform_tail_order.py` replays only the target-order arithmetic of
`DiagonalScale.exists_uniform_tail_order` for the Section 5 specialization.
The theorem bounds the scalar tail by

`2^(-J) q^(g(J+1)/2 - L)`

and with `g(j)=2*h*j` a requested target power `N` is reached once

`h*(J+1) - L >= N`.

The implementation selects the smallest such integer `J` above a requested
minimum order. `L` and `N` must be exact integers/Fractions; floating aliases are
rejected. The repository's runtime `h` follows the existing convention and is
interpreted exactly as its validated binary64 value through `Fraction.from_float`,
matching `slow_order_exact`. The certificate rechecks both the physical form and
the displayed formal gain condition `g(J+1) >= 2*(N+L)`.

## Truth boundary

This increment **does not** assert the hypotheses needed to apply the infinite
tail theorem to the actual reconstructed PDE background. In particular it does
not prove:

- existence/smoothness of hierarchy coefficients at every natural order;
- the stagewise actual `C[j,m]` rows for every order;
- an infinite coherent `exists_admissibleScales` schedule for the reconstructed hierarchy;
- exact retained recurrence cancellation at every order;
- that the actual PDE residual is bounded by the scalar dyadic tail;
- all-jets-flatness, Proposition 5.3, or super-algebraic convergence.

Accordingly the runtime certificate reports
`all_order_hierarchy_verified=false`,
`infinite_diagonal_schedule_verified=false`,
`pde_residual_tail_verified=false`, `all_jets_flat=false`,
`super_algebraic=false`, and `paper_exact=false`.

## Upstream handoff to Agent 2

The coefficient-level blocker is now stated in the pinned theorem's actual
quantifiers: Agent 2 must expose a hierarchy-owned construction **for every
requested coefficient order**, with smooth/analytic derivative and support data
sufficient to instantiate the stagewise template bounds, and exact retained
recurrence-row decomposition on the same hierarchy/source/coefficient state.
Agent 7 does not need a guessed global numeric supremum over all stages; it needs
an all-order provider whose totality and stagewise ownership are justified.

The currently landed finite eta/mixed-derivative frontier (including PR #301's
`partial_eta^3(Omega_k/X)`) is therefore genuine progress but cannot yet satisfy
that universal hierarchy quantifier.
