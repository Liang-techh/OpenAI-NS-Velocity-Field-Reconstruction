# Schedule-pressure scalar provenance

Pinned formal source: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

Relevant modules:

- `NavierStokes/FlatCutoff.lean`: `FlatCutoff.edge 1 x = 0` for `x <= 0`, otherwise `exp(-1/x^2)`.
- `NavierStokes/OutgoingSchedule.lean`: `sigma(x) = edge(1,x)/(edge(1,x)+edge(1,1-x))`.
- `NavierStokes/OutgoingTail.lean`: obtains an existential global derivative bound, stores it as noncomputable `stepBound := Classical.choose ...`, and sets `flattenLength = 10*(stepBound+1)*log 2 + 1`.
- `NavierStokes/SchedulePressure.lean`: `shapeExponent(d,y) = 1 - sigma((y-d.core.endpoint)/flattenLength)`.

## Constructive derivative-bound witness

For `0 < x < 1`, set

`a = exp(-1/x^2)` and `b = exp(-1/(1-x)^2)`.

Direct differentiation gives

`sigma'(x) = 2ab/(a+b)^2 * (x^-3 + (1-x)^-3)`.

The formula is symmetric under `x -> 1-x`, so it suffices to treat `0 < x <= 1/2`. On that interval, `b >= exp(-4)`, hence

`ab/(a+b)^2 <= a/b <= exp(4) a`.

Also `a <= exp(-4)` and, writing `t=1/x >= 2`, the function `t^3 exp(-t^2)` decreases for `t>=2`, so

`x^-3 exp(-1/x^2) <= 8 exp(-4)`.

Finally `(1-x)^-3 <= 8`. Therefore

`sigma'(x) <= 2 exp(4) * (8 exp(-4) + 8 exp(-4)) = 32`.

Outside `[0,1]`, `sigma` is locally constant, so the derivative is zero. Thus `S=32` is a valid explicit witness for the existential derivative bound used by `OutgoingTail.lean`, and this repository takes

`L_explicit = 10*(32+1)*log 2 + 1`.

## Truth boundary

This increment is an executable constructive replacement for an opaque theorem-side choice. It does **not** assert that Lean's particular `Classical.choose` term is definitionally equal to `32`, and it does not yet materialize the complete outgoing `finalAngular`, `clockWeight`, or `SchedulePressure.axisPressure`. Consequently Stage 1 remains `formal-structure`; `paper_exact_velocity_available` must remain false.
