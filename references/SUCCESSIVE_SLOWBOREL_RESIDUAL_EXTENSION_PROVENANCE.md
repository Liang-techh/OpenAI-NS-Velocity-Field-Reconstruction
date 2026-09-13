# Section 5 successive SlowBorel residual-extension provenance

## Scope

This increment couples three already fail-closed finite-prefix layers on Agent 7's convergence lane:

1. hierarchy/source-bound normalized-template constants `C[j,m]`;
2. the pinned recursive SlowBorel/DiagonalScale cutoff schedule; and
3. exact retained-recurrence cancellation followed by the first-omitted-order residual majorant.

`background_prefix_extension.py` first certifies one successive truncation step `N -> N+1`. The extension is accepted only when the new hierarchy-bound certificate literally preserves every old `C[j,m]` witness, the new schedule preserves every old local/recursive scale, and the new exact cancellation certificate literally preserves every old retained identity. Cutoff and recurrence evidence must use one hierarchy and source revision, and the recurrence prefix must stay on one coefficient-state identity.

The two residual majorants must be evaluated at the same `q`, `h`, and physical base power. The first omitted exponent is then checked with exact `Fraction` arithmetic and must gain exactly `2*h`. The code does **not** infer improvement from that exponent alone: the extended omitted coefficient prefactor may grow, so the new certified full-residual log-majorant must actually be no larger at the same `q`.

The same module now also composes any caller-supplied nonempty **finite** tuple of these one-step witnesses. Adjacent steps are accepted only when the earlier extended cutoff witness is exactly equal to the later previous cutoff witness and the earlier extended residual witness is exactly equal to the later previous residual witness. The whole tuple must therefore stay on one hierarchy/source/coefficient-state identity and one `q/h/base_power` regime, advance exactly one truncation order per step, gain exactly `2*h` in the first-omitted exponent per step, and retain a nonincreasing endpoint residual log-majorant. This blocks splicing two independently valid but incompatible finite-prefix paths.

## Mathematical boundary

This remains a finite-prefix implication only. It does not provide the missing paper hierarchy, does not derive uniform all-order `C[j,m]`, and does not prove that the finite chain can be extended for every `N`. In particular, accepting arbitrary *chosen finite lengths* is not an infinite coherent family and does not establish Proposition 5.3, all-jets-flatness, or super-algebraic decay.

In particular:

- fixtures cannot be promoted to paper data;
- numerical near-zero recurrence values never enter the cancellation gate;
- changing any previously certified bound, dependency, theorem identity, source revision, coefficient state, cutoff scale, `q`, `h`, or base exponent fails closed;
- an exponent gain with a larger finite residual majorant also fails closed;
- adjacent finite witnesses with different shared cutoff/residual endpoints fail closed;
- `paper_exact` and `infinite_coherent_family` are hard-coded `False` on the finite-chain witness.

## Verification

`tests/test_background_prefix_extension.py` checks a valid one-order extension, exact `2*h` first-omitted exponent gain, recursive-scale prefix preservation, strict finite-majorant improvement for a controlled fixture, and fail-closed rejection of an altered old `C[j,m]`, an altered old exact recurrence identity, a changed `q`, and a coefficient-prefactor increase large enough to defeat the slow-order gain.

It also checks a coherent two-step `1 -> 2 -> 3` chain, exact cumulative `4*h` first-omitted exponent gain, endpoint residual improvement, rejection of two individually valid steps whose shared residual endpoint was cross-wired through a different `q`, and rejection of an empty chain. The fixture uses exact rational formal cancellation atoms and explicitly supplied finite coefficient majorants. It is not paper coefficient data.

## Upstream/downstream handoff

Agent 2 must still emit genuine hierarchy-owned coefficient rows together with exact retained recurrence decompositions and the analytic derivative/support information needed to prove real `C[j,m]` bounds. More specifically, turning this finite-chain machinery into the manuscript's arbitrary-order statement requires a hierarchy-owned construction which, for every requested order `N`, supplies the next exact recurrence row and derivative/support bound family under one uniform analytic envelope. Agent 7 can then populate the existing provenance gates and prove a uniform extension theorem rather than merely compose already-supplied finite witnesses.

The next convergence milestone is therefore not another longer fixture chain. It is a constructive/uniform extension certificate over the real hierarchy that proves every finite prefix is available with the paper inequalities, from which the infinite all-jets-flat/super-algebraic closure can be derived. The finite coherent-chain witness added here deliberately stops before that theorem.
