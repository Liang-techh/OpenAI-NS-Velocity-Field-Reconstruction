# First-omitted physical tail power provenance

Issue: #2  
Lane: Agent 7 — downstream all-order analytic closure  
Formal source: `openai/NavierStokesAndEuler@f9e8bc5b38b6e212696e8a30e3e91517af887bbd`

## Purpose

This increment binds the finite retained-recurrence prefix already owned by the
hierarchy provider to the first omitted slow order and the exact ordinary /
physical q-power arithmetic used by the pinned SlowBorel tail theorems.  It is
stacked on the joint physical-prefix and exact uncut-radius chain; it does not
implement Agent 2 hierarchy jets or coefficient recurrences.

The formal source keeps three levels distinct:

1. `SlowBorelBase.exists_stage_blown_power_bound`: a fixed positive stage `j`
   retains the full slow power `q^(2*h*j)`.
2. `SlowBorelBase.ordinary_tail_bound`: after the recursive diagonal schedule
   absorbs coefficient/cutoff growth, the tail after prefix `J` has factor
   `(1/2)^J * q^(h*(J+1)-m)` for ordinary derivative order `m`.
3. `SlowBorelBase.exists_physical_uncut_tail`: physical-chart composition
   spends one further q-power per derivative; the theorem invokes the ordinary
   tail with target `P+m` and obtains a physical `q^P` bound after also supplying
   the independent finite chart constant.

The executable ledger therefore derives, from one selected prefix only,

- first omitted order `J+1`;
- full fixed-stage power `2*h*(J+1)`;
- diagonal tail power `h*(J+1)`;
- ordinary derivative power `h*(J+1)-m`;
- physical derivative power `h*(J+1)-2*m`;
- the same dyadic prefactor `2^-J` and exact uncut radius inherited from the
  provider-owned recursive schedule.

For the regression request `h=1/8`, `M=1`, `P=1/8`, `Jmin=4`, the selected
prefix is `J=16`, so the first omitted order is `17`.  The full stage power is
`17/4`; the diagonal tail power is `17/8`; the ordinary powers are `(17/8,9/8)`
and the physical powers are `(17/8,1/8)`.  The shared dyadic prefactor is
`2^-16=1/65536`.

## Exact cancellation boundary

The input dependency chain certifies exact retained recurrence identities only
for the retained prefix `0..J`; positive orders `1..J` are also tied to the same
own-order coefficient artifacts used by every `C[j,m]` majorant and support
witness.  This increment deliberately stops at `J`.  It does **not** synthesize
or assume a cancellation identity at the first omitted order `J+1`.

Consequently the nonzero/unknown defect boundary is explicit: the q-power of
the first omitted order is exposed, but its coefficient prefactor is not
claimed zero or bounded by this new certificate.

## Truth boundary

This increment verifies finite exact power bookkeeping only.  The following
remain false/open:

- first-omitted coefficient majorant verified;
- first-omitted recurrence cancelled;
- physical-chart finite constant verified;
- actual PDE residual evaluated or bounded;
- infinite `slowSum` materialized as an executable global object;
- total all-order hierarchy verified;
- infinite diagonal schedule verified;
- Proposition 5.3 / all-jets-flat closure;
- super-algebraic residual convergence;
- paper-exact reconstruction.

Passing regression or CI cannot promote any of these statements.

## Upstream requirement returned to Agent 2

The remaining coefficient-level blocker is still a total hierarchy-owned
assembled coefficient provider.  For every positive order demanded by an
arbitrary physical target request it must produce the actual assembled
coefficient object together with theorem-backed `C[j,m]` majorants, the common
compact-support witness, and exact retained-recurrence evidence under one
coefficient state.  The current finite PositiveAxis Picard/eta-jet frontier
cannot yet discharge that universal contract.
