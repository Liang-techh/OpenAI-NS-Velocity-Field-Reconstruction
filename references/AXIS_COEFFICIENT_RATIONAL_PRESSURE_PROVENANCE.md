# Exact rational reference source and pressure factors

Status: chosen-kernel coefficient arithmetic, not a paper-exact Stage 1 constructor.

## Source and scope

The pressure operator is pinned to
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`,
`NavierStokes/AxisContraction.lean`, `naturalRemainder`, with the coefficient
maps in `AxisOperators.lean` and `AxisWeightEstimates.lean`. The existing
`AXIS_COEFFICIENT_WIDE_NATURAL_PRESSURE_PROVENANCE.md` records this mapping.
No new Lean build or theorem verification is claimed here.

`axis_coefficient_rational_pressure.py` evaluates the same reference source
and linear pressure operator using exact `Fraction` arithmetic. It accepts
one actual amplitude state and reconstructs its reference angular coefficients
from `RationalAxisCoefficientData`; it does not reuse a possibly customized
angular-square provider from a manually constructed wide source state.

The selected binary64 h, j, sigma, and evaluation eta are interpreted as exact
binary rationals. A is computed as 1/2+h in rational arithmetic. Lambda is
the exact rational value of the selected finite Decimal. The coefficient
window is [-11/10,11/10], checked after conversion to binary64. This scope
does not certify earlier schedule quadratures or the selection of those inputs.

## Coefficient derivation

Radial indices are ordinary power-series coefficients. Parameter indices
are full derivatives, so only parameter products carry binomial factors.
Let phi[n,m] be the reference angular jet and
B[k] = (a^2)^(-1) d_eta^k(a^2). The rational amplitude Bell recurrence
computes B exactly from the normalized amplitude gradient.

The normalized source S = a^(-2) d_eta^m [a^2 phi^2]_n is

```
S[n,m] = sum(k=0..m) binom(m,k) B[k]
           * sum(i=0..n) sum(l=0..m-k)
               binom(m-k,l) phi[i,l] phi[n-i,m-k-l].
```

The pinned pressure is

```
j1(-4 A eta * primitive(source)
   + d * parameterPrimitive(source) - 2 eta * mulY(source)).
```

Primitive contributes S[n-1,m]/n; parameterPrimitive contributes
S[n-1,m+1]/n, including differentiation of a^2; mulY contributes
S[n-1,m]. All three vanish in radial row 0. The final j1 takes the
predecessor input row and divides by n^2, with row 0 equal to zero.

An independent simplification uses the fact that eta and d have only radial
row 0. Thus P[0,m]=P[1,m]=0 and, for n>=2,

```
P[n,m] = sum(k=0..m) binom(m,k)
           * (d^(k) S[n-2,m-k+1]
              - (4 A + 2(n-1)) eta^(k) S[n-2,m-k])
         / ((n-1)n^2).
```

In particular phi[0]=1 gives S[0,m]=B[m] and
P[2,0]=(d B[1]-(4 A+2)eta)/4. These identities expose both radial
index shifts and the derivative of the amplitude, independently of the
implementation's composed operator calls.

The normalized first parameter derivative is
P[2,1]=(d' B[1]+d B[2]-(4 A+2)(1+eta B[1]))/4.
The eta B[1] term is necessary: these are derivatives of the full
a^2-scaled pressure, followed by division by a^2, rather than derivatives
of the normalized factor alone.

## Recorded verification

The focused pressure test file passed: **4 passed in 0.76 seconds**;
compilation also passed. Coverage includes direct polynomial H/H' and
Bell B[1]/B[2] identities at eta=0 and eta=0.01, the first two nonconstant
radial source coefficients, pressure radial shifts and parameter product
rules, signed log bounds, and invalid controls before zero shortcuts.
The comparison with the existing Decimal evaluator is diagnostic only.
No full-suite result is attributed to this increment.

## Log enclosure and remaining boundary

For a nonzero exact factor F, `rational_log_enclosure(abs(F))` supplies
[l,u]. If the common amplitude source has log bounds [L,U], the returned
signed coefficient has exact sign(F) and log magnitude in [2L+l,2U+u].
Its width is 2(U-L)+(u-l). The common amplitude is retained once, so
correlations are preserved through the complete finite factor calculation.
Zero factors have sign 0 and no finite logarithm. The public result
dataclass is a container conditional on its supplied factor and interval;
the actual-state factory computes them through the rational path.

This removes float/Decimal rounding from these finite normalized source
and pressure factors. It does not enclose subsequent mixed-scale sums,
the axial zStar dependency, all Picard iterates, or a global weighted
coefficient norm. Stage 1 and the full velocity construction remain
incomplete. No metadata or passing test promotes their status.
