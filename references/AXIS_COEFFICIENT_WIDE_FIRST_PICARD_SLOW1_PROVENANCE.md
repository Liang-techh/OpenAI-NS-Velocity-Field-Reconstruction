# Genuine first-Picard `slow1(x1)` provenance

This increment is pinned to `openai/NavierStokesAndEuler` commit
`f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.

The formal source defines `bu = O.average u`,
`averageCoefficient(d) = (2*d.D) • d.eta`,
`angularSlowCoefficient(d) = (2*d.h) • d.eta`, and

`slow1 = j2(product(product(averageCoefficient d,bu) + product(angularSlowCoefficient d,u),phi)) + param2 bu (product d.d phi) + dot2 (product (averageCoefficient d) bu) phi + product d.d (mixed2 bu phi) - param2 phi (product d.d u)`.

`AxisOperators.averageData` gives the exact coefficient row factor `1/(n+1)`. `mixed2 = inverseMixed 2`; pinned `inverseMixedJet` places one eta derivative on its first argument, the radial-Euler factor of the second argument, and the regular-inverse divisor on the output row. Together with `j2`, `dot2`, and `param2`, each nonzero output row `n` uses `radialDivisor(2,n-1)=n(n+1)`. Row zero is therefore the exact zero-datum inverse row.

The implementation consumes only the landed `ActualScheduleWideFirstPicardState` and its actual SchedulePressure AxisData. It keeps ordinary `Lambda^0..Lambda^-4` and pressure-linear `a(eta)^2 Lambda^-1..Lambda^-3` numerators separate. Pressure derivatives use the landed normalized pressure jet chain; the common amplitude is not narrowed to binary64.

## Truth boundary

This certifies only executable coefficient-level `slow1(x1)`. It does not certify global weighted `AxisCoefficientSpace` membership, complete `naturalRemainder(x1)`, `x2`, Picard convergence, a fixed point, `NaturalProfileAssembly`, or paper-exact velocity. Structural regression tests do not promote any of those claims.
