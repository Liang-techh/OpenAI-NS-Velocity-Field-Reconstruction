# Section 10 endpoint-majorant ladder provenance

Status: **formal-structure** only. This file does not certify a paper-exact force or velocity.

## Purpose

`src/openai_ns_reconstruction/section10_endpoint_ladder.py` packages the already-landed one-degree endpoint-limit implication into a contiguous finite family of derivative degrees. The purpose is to make the Section 10 smooth-extension prerequisite fail closed at the interface level before the actual residual supplies all-order estimates.

The pinned construction requires locally uniform limits of every spacetime derivative of the closed-past Navier--Stokes residual as `t -> 1-`. The repository currently represents a sufficient hypothesis for one derivative degree by `EndpointPowerLawMajorant`: an independently proved bound

`sup_{x in K_m} ||partial_t D^n R(t,x)|| <= C_n (1-t)^(-alpha_n)`, with `alpha_n < 1`.

For a finite ladder `n=0,...,N`, the new adapter requires every degree exactly once, one common spatial window, and passes every member through the existing `section10_endpoint_localization_transfer`. Therefore every recorded bound must have endpoint `T=1` and begin in the official `t>=3/4` unit plateau, where `chi=1`, `chi'=0`, and `chi^2-chi=0` exactly for the executable Section 10 switch representative.

On the common validity interval, the adapter computes the maximum of the exact antiderivative Cauchy/tail budgets. Conditional on the supplied hypotheses, that gives one finite-order modulus controlling derivative degrees `0..N` simultaneously.

## Independent regression

`tests/test_section10_endpoint_ladder.py` does not estimate the majorants from residual samples. It independently evaluates the closed-form antiderivative

`C/(1-alpha) * (1-t)^(1-alpha)`

and the corresponding two-time integral for each supplied degree, then compares those values with the ladder output. Additional tests verify fail-closed behavior for missing/duplicate derivative degrees, mixed spatial windows, transition-collar start times, and non-`t=1` endpoints.

## What this closes

- one finite family can no longer silently omit an intermediate derivative degree;
- one common endpoint window is derived from the latest member validity start;
- every degree is gated through the official late time-localization plateau;
- finite-order endpoint Cauchy/tail budgets can be queried with one common modulus.

## What remains unresolved

This is **not** evidence that the actual Section 9/10 residual satisfies any of the supplied power-law majorants. It does not construct endpoint derivative limits, prove the full-spacetime tensor compatibility required by `CandidateFromLimits`, prove all orders at once, choose the Borel scale from genuine analytic template bounds, prove past/future jet matching, or certify a smooth compactly supported force through `t=1`.

The paper-exact gate must remain false until those upstream hypotheses are derived from the completed local field and the infinite family of endpoint derivatives is closed.
