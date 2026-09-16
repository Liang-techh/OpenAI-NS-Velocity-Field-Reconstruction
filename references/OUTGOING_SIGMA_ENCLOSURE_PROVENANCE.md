# Rational enclosures for the actual outgoing smooth step

Status: arithmetic prerequisite for certified schedule pressure integration.
This does not certify the complete pressure integral or Stage 1.

## Pinned function and the missing dependency

The source is `OutgoingSchedule.sigma` and `FlatCutoff.edge` at
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.
The repository's `schedule_pressure.py` records the defining formula:

```
edge(x) = exp(-1/x^2) for x>0, and 0 otherwise
sigma(x) = edge(x)/(edge(x)+edge(1-x)).
```

This differs from `cutoffs.smooth_step`, whose exponent has first powers.
The existing `outgoing_tail._sigma_primitive` uses Gauss quadrature inside
the schedule amplitude and clock weight. Thus bounding only the outer
pressure quadrature would leave an uncontrolled inner integral.

`outgoing_sigma_enclosure.py` supplies rational bounds for the actual
sigma and its primitive. Inputs are exact `Fraction` values, not silently
rounded floating-point nodes. No alternative bump or fitted polynomial is
introduced. The pinned opaque step-bound choice is not changed.

## Exponential bound

For rational x>=0, choose s>=0 such that y=x/2^s<=1. The positive Taylor
sum T_N of exp(y), through degree N, has remainder bounded above by

```
R_N = (y^(N+1)/(N+1)!)/(1-y/(N+2)).
```

The ratio of each successive omitted term is at most y/(N+2)<1. Hence
[1/(T_N+R_N),1/T_N] encloses exp(-y). Let eps=min(requested_width,1)
and refine until R_N<=eps/(4*2^s). Choose a dyadic grid step
delta<=eps/(8*2^s), round interval endpoints outward, and square the
interval s times, rounding outward after each operation. Intersecting with
[0,1] is valid because exp(-x) lies there.

On [0,1], squaring multiplies interval width by at most 2, and outward
rounding adds at most 2delta. Including initial rounding, the final width
is bounded by

```
2^s R_N + 2delta*2^s + 2delta*(2^s-1) <= 3eps/4.
```

All operations and comparisons use rational arithmetic. Caps on Taylor
terms and squarings must fail explicitly when sufficient refinement is not
obtained; a cap is not an error estimate.

## Stable sigma and convergent primitive

For 0<x<=1/2, let t=1/x^2-1/(1-x)^2>=0. Then
sigma(x)=r/(1+r), where r=exp(-t). This increasing map has derivative
at most 1 for r>=0, so it preserves the requested width bound.
For x>1/2, use sigma(x)=1-sigma(1-x). The values outside (0,1) and
at 1/2 are exact.

Inside (0,1), the derivative is

```
sigma'(x)=2 sigma(x)(1-sigma(x)) (x^-3+(1-x)^-3)>0.
```

Consequently, left and right endpoint sums provide rigorous integral
bounds without a sampled derivative maximum or an empirical quadrature
error. The primitive I(x)=integral_0^x sigma obeys I(x)=0 for x<=0,
I(x)=x-1/2 for x>=1, and I(x)=x-1/2+I(1-x) for x>1/2.
Only the interval [0,y], with y<=1/2, needs refinement.

For N equal cells, bound each cell by its width times sigma's enclosed
left and right endpoint values. Endpoint intervals are rounded outward
onto one common dyadic grid to control denominator growth. With each
sigma interval width at most eps/(16y) before rounding and grid step at
most eps/(16y), the total endpoint uncertainty contributes at most
3eps/16 to the integral width. The remaining Darboux gap tends to zero
as N increases (it is bounded by y/N). Thus the scheme converges below
the requested width when the cell cap permits it. Returned bounds are
accepted only after their exact width passes the requested tolerance.

## Remaining pressure work

These bounds close the scalar smooth-step integral dependency. The actual
clock-weight formula still needs interval evaluation, including its derived
transition geometry and constants. Then the outer pressure integral and
its parameter derivatives need certified refinement before zStar and the
axial reference coefficients can consume pressure bounds. The existing
floating-point pressure path is not automatically certified by this module.

## Recorded verification

The focused test file passed **5 tests in 0.21 seconds**; both new Python
files compiled. Tests cover exact branches, rational interval widths and
caps, high-precision exponential/sigma diagnostics, sigma and primitive
reflection, and primitive refinement at x=1/2 with widths 1/64 and 1/256.
The existing Gauss-rule comparison is diagnostic, not the proof of the
enclosure. No full-suite or complete pressure certification is claimed.

## Next integral reductions

Two further reductions avoid differentiating sigma in those future bounds.
Write k=1-h>0 and rho for the actual tail parameter. Integration by parts
gives the exact tail debt identity

```
tailDebt = rho/(1-rho)
           * (exp(3k)-2k*integral_0^1 exp(k*(1+2x))*sigma(x) dx).
```

The integrand is increasing and positive, so the same endpoint-enclosure
strategy applies. Cancellation must be retained as interval subtraction;
positivity cannot be asserted unless the resulting lower bound proves it.

For release lag, put A(t)=release_rate_primitive(t), T=ramp_end, and
q0=(lambda-h)/(1-lambda). Since release_source=(1-h)-A', integration
by parts gives

```
releaseLag(T) = exp(-A(T))
                * (q0+1+(1-h)*integral_0^T exp(A(t)) dt) - 1.
A(T) = 1-(lambda+h)/2.
```

Here A'>=0 follows from release_slope>=-1, so exp(A) is increasing.
On the long middle interval A is exactly (1-lambda)/2; its integral
is that interval's length times exp((1-lambda)/2). Only the two unit
transitions remain to be integrated. These are mathematical reductions
for the next implementation; the current module does not compute either
tailDebt or releaseLag.
