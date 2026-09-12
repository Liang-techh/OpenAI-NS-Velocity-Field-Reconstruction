# Section 5 first-Picard second-parameter-jet provenance

## Source

Primary paper: *Finite Time Blowup for Navier–Stokes* (OpenAI, September 2026).

Paper PDF: `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`

Relevant printed pages: 45-52, especially Lemma 5.1 / Eq. (5.7).

For

`(G g)_i(xi,eta)=integral_0^xi (s/xi)^c_i g_i(s,eta) ds`,

with `c=(0,0,2,0,3,1)`, the Volterra kernel has no eta dependence. Under the smoothness used in Lemma 5.1, differentiating under the integral twice gives

`partial_eta^2(G g)=G(partial_eta^2 g)`.

Applied to `W_n^(0)=G f_n`, this is the exact second-parameter commutation needed before differentiating the genuine `k=1` Picard right-hand side, whose `A1 partial_eta W_n^(0)` row introduces `partial_eta^2 W_n^(0)`.

## Executable mapping

`src/openai_ns_reconstruction/background_first_picard_second_parameter_jet.py` adds `first_picard_second_parameter_jet_eq_5_7(...)`. It reuses the landed first-Picard value/first-parameter path for `(G f_n, partial_eta G f_n)` and applies the same Eq. (5.7) singular inverse directly to an explicit analytic `partial_eta^2 f_n` provider for the new row. No production finite difference, sampled fit, generic cutoff, or independently supplied Picard iterate is introduced.

`tests/test_background_first_picard_second_parameter_jet.py` uses an independent polynomial-in-radius, cubic-in-eta source. It checks all three rows against closed-form singular integrals and checks the second derivative separately against a centered eta difference of the previously landed first-parameter value path. The finite difference is test-only and the fixture is not manuscript coefficient data.

## Truth boundary

This increment is **formal-structure / solver infrastructure**. The new `partial_eta^2 f_n` provider is explicit caller data. A sampled, fitted, or otherwise caller-authored second derivative is not paper-exact provenance.

The current strong repaired-history path owns `f_n` and `partial_eta f_n`, but it does not yet own the additional mixed `phi/U/beta` jets required to differentiate the complete lower-history forcing a second time. In particular this increment does not claim hierarchy-owned `partial_eta^2 f_n`, hierarchy-owned `partial_eta^2 W_n^(0)`, `partial_eta W_n^(1)`, another Picard iterate, convergence of the Lemma-5.1 series, or a recursively materialized final coefficient.

Stage 2 remains `formal-structure`; `paper_exact_velocity_available=false`. The real leading mixed profile remains upstream in Issue #1, and common-strip/(5.17) bounds, completion of the infinite recursive cutoff argument, and Proposition 5.3 all-jets residual decay remain open.
