# Section 5 successive SlowBorel residual-extension provenance

## Scope

This increment couples three already fail-closed finite-prefix layers on Agent 7's convergence lane:

1. hierarchy/source-bound normalized-template constants `C[j,m]`;
2. the pinned recursive SlowBorel/DiagonalScale cutoff schedule; and
3. exact retained-recurrence cancellation followed by the first-omitted-order residual majorant.

`background_prefix_extension.py` certifies one successive truncation step `N -> N+1`. The extension is accepted only when the new hierarchy-bound certificate literally preserves every old `C[j,m]` witness, the new schedule preserves every old local/recursive scale, and the new exact cancellation certificate literally preserves every old retained identity. Cutoff and recurrence evidence must use one hierarchy and source revision, and the recurrence prefix must stay on one coefficient-state identity.

The two residual majorants must be evaluated at the same `q`, `h`, and physical base power. The first omitted exponent is then checked with exact `Fraction` arithmetic and must gain exactly `2*h`. The code does **not** infer improvement from that exponent alone: the extended omitted coefficient prefactor may grow, so the new certified full-residual log-majorant must actually be no larger at the same `q`.

## Mathematical boundary

This is a finite-prefix implication only. It does not provide the missing paper hierarchy, does not derive uniform all-order `C[j,m]`, does not prove the existence of an infinite coherent chain of certificates, and does not establish Proposition 5.3, all-jets-flatness, or super-algebraic decay.

In particular:

- fixtures cannot be promoted to paper data;
- numerical near-zero recurrence values never enter the cancellation gate;
- changing any previously certified bound, dependency, theorem identity, source revision, coefficient state, cutoff scale, `q`, `h`, or base exponent fails closed;
- an exponent gain with a larger finite residual majorant also fails closed;
- `paper_exact` is hard-coded `False`.

## Verification

`tests/test_background_prefix_extension.py` checks a valid one-order extension, exact `2*h` first-omitted exponent gain, recursive-scale prefix preservation, strict finite-majorant improvement for a controlled fixture, and fail-closed rejection of an altered old `C[j,m]`, an altered old exact recurrence identity, a changed `q`, and a coefficient-prefactor increase large enough to defeat the slow-order gain.

The fixture uses exact rational formal cancellation atoms and explicitly supplied finite coefficient majorants. It is not paper coefficient data.

## Upstream/downstream handoff

Agent 2 must still emit genuine hierarchy-owned coefficient rows together with exact retained recurrence decompositions and the analytic derivative/support information needed to prove real `C[j,m]` bounds. Agent 7 can then populate the existing provenance gates and this successive-extension certificate with real hierarchy data.

The next convergence milestone is an arbitrary-length coherent extension family whose bounds are generated from the actual hierarchy and whose recursive schedule/tail estimates imply the manuscript's infinite all-jets-flat/super-algebraic closure. Finite successive witnesses alone are not that theorem.
