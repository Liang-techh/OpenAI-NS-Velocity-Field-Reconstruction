# Section 5 first-Picard radial-jet provenance

## Source

Primary paper: *Finite Time Blowup for Navier–Stokes* (OpenAI, September 2026).

Paper PDF: `https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf`

Relevant printed pages: 45-52, especially Lemma 5.1 / Eq. (5.7).

The landed Eq. (5.7) solver uses

`W_n=(phi_n,U_n,K_n,Pi_n,d_xi phi_n,d_xi U_n)`

and the singular inverse

`(G g)_i(xi)=integral_0^xi (s/xi)^c_i g_i(s) ds`,

with `c=(0,0,2,0,3,1)`.  For the first Picard term `W^(0)_n=G f_n`, direct differentiation at positive radius gives

`d_xi W^(0)_i = f_{n,i}(xi,eta) - c_i W^(0)_i/xi`.

This identity is a consequence of the displayed singular inverse and does not introduce a fitted derivative or a new manuscript assumption.

## Executable mapping

`src/openai_ns_reconstruction/background_first_picard_radial_jet.py` adds `first_picard_radial_jet_from_hierarchy(...)`.

The function accepts only a `Section5LowerHistoryJetHierarchy`, so `f_n` is generated through the already-landed contiguous hierarchy, analytic Eq. (5.2) beta path, `actualLowerSource`, and exact displayed Eq. (5.7) forcing.  It evaluates the existing genuine first term `G f_n`, then obtains its positive-radius radial derivative from the differentiated singular-integral identity above.  It also records the pointwise singular-ODE defect

`d_xi W + xi^-1 diag(c) W - f_n`.

The adapter deliberately requires `xi>0`.  It does not infer an axis derivative from a single source value, because doing so would require a neighborhood regularity statement.  It also does not manufacture `d_eta W`; therefore it does not yet unlock the full next Picard application `G(A0 W + A1 d_eta W + f_n)`.

`tests/test_background_first_picard_radial_jet.py` independently differentiates the already-landed integral term with a centered radial finite difference and compares that result with the analytic derivative returned by the new adapter.  A separate check evaluates the singular ODE residual using a fresh hierarchy-generated forcing callback.  The test fixture is analytic and is only a regression oracle; it is not promoted to a paper coefficient.

## Truth boundary

This increment is **formal-structure / solver infrastructure**, not a completed positive-order coefficient and not a paper-exact background.

The first Picard term is only the `k=0` term of the Lemma-5.1 series.  The coefficient matrices start contributing at the next Picard application, which still requires a controlled eta derivative of the current iterate.  The order-zero leading profile remains upstream from Issue #1, and positive-order hierarchy data have not been recursively materialized to all orders.

No common analytic strip, coefficient bounds (5.17), full Picard convergence certificate, recursive cutoff local-finiteness proof, Proposition 5.3 all-jets residual decay, or final cutoff-summed background is established here.  `paper_exact_velocity_available` remains false.
