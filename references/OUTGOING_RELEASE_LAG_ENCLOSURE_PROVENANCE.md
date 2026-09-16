# Outgoing release-lag enclosure

Status: certified scalar arithmetic prerequisite, not complete pressure or
Stage 1 reconstruction. Source: `OutgoingTail.releaseLag` at
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`;
the existing numerical counterpart is `TailData.release_lag_at_ramp_end`.

## Exact reduction

Interpret the actual TailData h and lambda as their exact binary rationals,
with 0<2h<lambda<1/10. Put c=1-lambda, k=1-h, a=c/2, B=(c+k)/2,
q0=(lambda-h)/c, and L=4 log(1/h). Let I(x)=integral_0^x sigma.
The ramp ends at T=L+2. Its rate primitive A satisfies A(T)=B.
Since release_source=k-A', integration by parts gives

```
lag(T) = exp(-B)*(q0+1+k*integral_0^T exp(A(t)) dt)-1.
```

The first unit transition has A(t)=c(t-I(t)). The middle interval has
length L and constant A=a. On the final unit transition, with local
coordinate x, A=a+k I(x). Symmetry yields I(1-x)=I(x)-x+1/2,
so the first transition integral is exp(a)*integral_0^1 exp(-c I(x)) dx.
Consequently

```
F = integral_0^1 f(I(x)) dx,
f(v) = exp(-c v)+exp(k v),
lag(T) = exp(-B)*(q0+1+k*exp(a)*(L+F))-1.
```

This eliminates both explicit sigma derivatives and integration across
the potentially long middle plateau. No numerically evaluated geometry
or lag is used as a certified input.

## Shared-grid enclosure

The actual sigma is increasing, takes values in [0,1], and has integral
1/2 on [0,1]. On a uniform grid, outward-rounded sigma intervals give
left/right Darboux bounds for every prefix I(x_i) in a single cumulative
pass. Intersect with [0,1/2], with exact values at endpoints 0 and 1.
This replaces repeated independent integrations at each outer node.

For v>=0, f'(v)=k exp(kv)-c exp(-cv)>=k-c>0. Since I is increasing,
f(I(x)) is increasing. Evaluate f at the lower and upper prefix bounds
using rigorous exponentials, then form a second pair of endpoint sums.
Common dyadic rounding controls denominator growth without reversing any
inequality. Prefix and outer Darboux errors both vanish under refinement;
arithmetic precision is budgeted separately.

The logarithm in L uses `rational_log_enclosure(1/h)`. Exponentials use
the validated rational negative exponential and the bounded positive
exponential helper. Products are positive until the final subtraction of
1, which shifts both endpoints exactly. Accept only a positive final lag
interval whose exact width meets the requested tolerance. Reaching a cap
without those conditions is explicit failure, not an approximate answer.

## Recorded verification

The focused file passed **3 tests in 0.37 seconds**; compilation passed.
The actual fixture P=2, m=1, lambda=0.05, wait=30, h=0.01 used 128 cells
for requested width 1/256, returning width about 0.002458971454861.
Requested width 1/64 used 32 cells. Both positive intervals overlapped
and enclosed the legacy numerical lag as a diagnostic comparison.
Insufficient cell caps and invalid controls failed explicitly. These
focused results do not establish full-suite or complete pressure acceptance.

## Remaining scope

The complete pressure calculation still requires combining this lag with
tailDebt to enclose decayHold, constructing certified transition geometry
and clock weights, and integrating the pressure kernel and its parameter
derivatives. The constructor does not silently certify those consumers.

For the next geometry step, positive lag bounds [qL,qU] and tail-debt
bounds [dL,dU], with qL>dU, give

```
decayHold in [log(qL/dU)/(1-h), log(qU/dL)/(1-h)].
```

Each endpoint logarithm must itself be enclosed with the rational logarithm
routine, taking its lower or upper endpoint respectively. Checking qL>dU
proves the required positive hold length; comparing interval midpoints would
not prove it. Both source constructors must use the same actual TailData h.
This displayed formula is a dependency specification, not an implemented
geometry certificate in the release-lag module.
