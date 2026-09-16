# Certified outgoing tail debt

Status: scalar parameter enclosure for the actual selected outgoing schedule.
The full pressure integral, axial reference, and Stage 1 remain incomplete.

## Pinned formula and inputs

The source is `OutgoingTail.tailDebt` and the outgoing sigma definition at
`openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`.
The executable point formula is in `outgoing_tail.TailData.tail_debt`.
For the existing constructive step bound S=32, define

```
rho = h exp(-5)/(16(32+1)) = h exp(-5)/528,
k = 1-h,
tailShape'(t) = (rho/2) sigma'((t-1)/2).
tailDebt = integral_1^3 exp(k t) tailShape'(t) dt / (1-rho).
```

The new constructor receives actual TailData and interprets its original
binary64 h as an exact rational. It recomputes rho through the certified
exponential; the rounded `data.rho` and numerical `data.tail_debt` are not
inputs to its certified arithmetic. The admissible range 0<h<1/20 implies
0<k<1 and 0<rho<1/10560. Other tail data affect later schedule geometry,
but are absent from this particular formula.

## Derivation

Substitute t=1+2x. The debt is rho/(1-rho) times
the Stieltjes integral of exp(k(1+2x)) against sigma on [0,1]. Integration
by parts, using sigma(0)=0 and sigma(1)=1, gives

```
J = integral_0^1 exp(k(1+2x))*sigma(x) dx,
F = exp(3k)-2k J,
tailDebt = rho/(1-rho) * F.
```

Both factors of J's integrand are nonnegative and increasing, so uniform
left/right endpoint sums enclose J. The sigma evaluator uses the pinned
squared exponent, as proved in `OUTGOING_SIGMA_ENCLOSURE_PROVENANCE.md`.
Endpoint interval products are rounded outward to a common dyadic grid
before summation. No floating-point quadrature error estimate is used.

For rational 0<=z<=3, e^z<27 follows from e<3. Enclose exp(-z), intersect
with [1/27,1], and invert endpoints. The inverse map has derivative magnitude
at most 729 there; an input width at most tolerance/729 therefore suffices
for the positive exponential's requested width.

The subtraction defining F can cause interval overestimation. Independently,
sigma' is nonnegative and has integral 1, so its Stieltjes integral is an
average of values in [exp(k),exp(3k)]. Consequently intersecting F with
[1,27] is valid and supplies a positive lower bound. This is a proved
intersection, not replacement of a missing coefficient or empirical clipping.

The map rho -> rho/(1-rho) is increasing on this domain. Multiply its
positive interval with the positive F interval. Return only when the exact
result width is at most the requested tolerance; reaching a cell or arithmetic
cap must fail explicitly. Doubling cell counts decreases the Darboux gap.

The refinement target is the final debt width, not J's width alone. Its
uncertainty is attenuated by 2k rho/(1-rho); imposing the final tolerance
directly on J would waste this small factor and may exhaust the cell cap.
Node arithmetic uses eta=min(1/1024,requested_width/1024). Each positive
exponential/sigma product has width at most 28eta before grid rounding;
its rounding adds at most eta/2. Meanwhile rho/(1-rho)<1/10000 and the
rho ratio uncertainty is at most eta/5280. Thus the residual arithmetic
uncertainty stays well below the final tolerance as cells are refined.

## Recorded verification

The focused file passed **3 tests in 0.68 seconds**; compilation passed.
For P=2, m=1, lambda=0.05, wait=30, h=0.01 (binary64 inputs interpreted
as specified above), the requested width 1/10^8 took 512 cells. The
returned rational interval, displayed approximately, was
[9.244482075814e-7, 9.340675096908e-7], with width about
9.619302109324e-9. Its exact rational width satisfied the tolerance.
The existing floating tail-debt value lay inside as a diagnostic check.
Tests also checked refinement overlap, the independent positive factor
bounds, and explicit rejection at an insufficient cell cap.
No full-suite result is attributed to this increment.

## Scope

This is a convergent arithmetic enclosure of one required actual parameter.
It does not replace the full tail's release-lag calculation, prove the
decay-hold inequality, certify transition geometry, or enclose the outer
pressure integral and its parameter derivatives. Those remain subsequent
dependencies. No new Lean build or full reconstruction claim is made.
