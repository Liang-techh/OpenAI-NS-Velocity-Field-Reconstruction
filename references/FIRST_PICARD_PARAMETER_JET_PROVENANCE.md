# Section 5 first-Picard parameter-jet provenance

## Source

Primary paper: *Finite Time Blowup for Navier–Stokes* (OpenAI, September 2026).

Paper PDF: `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`

Relevant printed pages: 45-52, especially Lemma 5.1 / Eq. (5.7).

For

`(G g)_i(xi,eta)=integral_0^xi (s/xi)^c_i g_i(s,eta) ds`,

with `c=(0,0,2,0,3,1)`, the kernel has no eta dependence.  Under the smoothness already assumed by Lemma 5.1, differentiation under the integral gives

`partial_eta(G g)=G(partial_eta g)`.

Applied to the first Picard term `W_n^(0)=G f_n`, this yields an exact-structure parameter derivative once a controlled analytic `partial_eta f_n` is available.

## Executable mapping

`src/openai_ns_reconstruction/background_first_picard_parameter_jet.py` adds `first_picard_parameter_jet_eq_5_7(...)`.  The value path calls the existing `first_picard_term_eq_5_7`; the parameter path applies the same landed singular inverse directly to an explicit `lower_order_source_parameter` callback.  Production does not finite-difference neighboring eta samples.

The interface deliberately does **not** call itself hierarchy-owned.  The current `Section5LowerHistoryJetHierarchy` supplies enough derivatives for the pointwise source `f_n`, including a second jet of `beta=V/X`, but differentiating the Eq. (5.6) `Omega/X` row once more requires higher regular-flux jets.  If those data are absent, this module requires the derivative callback rather than silently fitting one.

`tests/test_background_first_picard_parameter_jet.py` uses an analytic polynomial-in-radius, quadratic-in-eta source family.  It compares production against a separately written closed-form singular integral and also against an independent centered finite difference of the already-landed `G f_n` value path.  The fixture is a regression oracle only and is not manuscript coefficient data.

## Truth boundary

This increment is **formal-structure / solver infrastructure**.  A caller-supplied `partial_eta f_n`, even if numerically accurate, is not paper-exact provenance.  The increment does not claim that the repaired hierarchy yet owns the extra derivatives needed to construct that callback, and it does not execute the second Picard application with manuscript-derived data.

The next strict step is to lift the repaired lower-history data to the additional eta derivatives required by `partial_eta(Omega/X)` (or an equivalent analytically certified source-parameter provider), then compose this commutation with the landed `A0/A1` fields to form the genuine `k=1` Picard term.

No full Picard convergence certificate, recursively materialized all-order coefficient sequence, common analytic strip, coefficient bounds (5.17), infinite cutoff schedule, Proposition 5.3 all-jets residual decay, or final cutoff-summed background is established here.  `paper_exact_velocity_available` remains false.
