# Stage-1 complete mixed-scale `naturalRemainder(x0)` pair provenance

Status: **formal-structure**. This increment assembles the two already-landed
mixed-scale remainder branches at the actual contraction centre. It does not
claim a paper-exact leading profile, a Picard iterate, or a certified fixed
point.

## Pinned source

Mapped to `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`, principally:

- `NavierStokes/AxisContraction.lean`
  - `naturalRemainder`
  - `referencePair`
  - the angular `lin1`, `quad1`, `slow1` branch
  - the axial `lin2`, `slow2`, pressure branch
- `NavierStokes/AxisOperators.lean`
- `NavierStokes/AxisResolvent.lean`
- `NavierStokes/NaturalAxisCoefficients.lean`

At `x0 = referencePair` and theorem-selected `t = 1/Lambda`, the two output
components are the previously landed mixed-scale expressions

- angular: `naturalResolvent(inverseL * (lin1 + quad1 - slow1/Lambda))`;
- axial: `inverseL * (lin2 - slow2/Lambda + pressure)`, with the genuine
  `a^2 * phi0^2` pressure source propagated in signed-log form.

## Landed capability

`src/openai_ns_reconstruction/axis_coefficient_wide_natural_remainder.py`
introduces one typed value for the **complete pair**
`naturalRemainder(referencePair)`. The production factory accepts only the
actual `TailData`, `j`, and the existing numerical phase-resolution control,
then calls the landed theorem-grounded angular and axial builders.

Before exposing the pair, the type fails closed unless both halves agree on:

1. the exact schedule `TailData`;
2. the same schedule parameter `j`;
3. the same certified `sigma`;
4. the same theorem-selected coefficient-space `epsilon`;
5. the same positive Decimal `Lambda`.

No scale is collapsed during assembly. Angular `1/Lambda` remains an explicit
Decimal numerator/denominator decomposition. Axial `1/Lambda` remains the same,
and the pressure contribution retains its signed-log `a^2` scale. The pair's
`jet_pair(n,m,eta)` simply exposes the two typed decompositions at the same
coefficient coordinate.

There is no caller override for pressure datum, sigma, epsilon, `Lambda`, `C`,
amplitude, coefficient tables, resolvent cutoff, or derivative tables.

## Verification

The regression suite checks that:

- the assembled angular and axial halves share the same actual schedule datum,
  `j`, certified sigma/epsilon, and exact theorem-selected Decimal `Lambda`;
- `jet_pair` preserves the landed mixed-scale angular and axial jet objects,
  including the noncollapsed `1/Lambda` terms and signed-log axial pressure;
- pairing an axial state from a different **actual theorem-derived** schedule
  branch is rejected rather than silently mixed;
- non-`TailData` production inputs fail closed;
- provenance remains machine-readable with `formal-structure` and
  `paper_exact_velocity_available=false`.

The different-`j` axial state used by the compatibility regression is itself
built through `actual_schedule_reference_wide_axial_remainder_state`; it is not
a production surrogate or a replacement manuscript parameter.

## Explicit non-claims / next blocker

This increment completes the typed mixed-scale representation of
`naturalRemainder(x0)` only. It deliberately does **not** apply the fixed-point
map's outer `1/(2*Lambda)` factor. Therefore `Picard x1` is not yet
materialized.

The next smallest blocker is to carry that outer factor through both mixed
scales without binary64 collapse, forming the first genuine Picard iterate
`x1 = x0 + naturalRemainder(x0)/(2*Lambda)`, and only then advance toward
iterative fixed-point `phi/u`, derived average/pressure,
`NaturalProfileAssembly`, and support/moment/matching/cone validation.

No global all-index weighted `AxisSpace` membership/norm certificate or
paper-exact velocity field is claimed. `paper_exact_velocity_available=false`
remains mandatory.
