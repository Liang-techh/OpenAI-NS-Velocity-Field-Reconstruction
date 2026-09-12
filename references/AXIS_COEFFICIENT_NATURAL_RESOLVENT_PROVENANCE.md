# Axis coefficient natural resolvent provenance

## Scope

This increment materializes the pinned `AxisResolvent.naturalResolvent` action at the landed coefficient-jet representation level for the actual SchedulePressure construction. It does **not** upgrade Stage 1 to paper-exact and does not claim the global weighted `AxisSpace` norm/membership object.

## Pinned source

- Repository: `openai/NavierStokesAndEuler`
- Commit: `f9e8bc5b38b6e212696e8a30e3e91517af887bbd`
- Module: `NavierStokes/AxisResolvent.lean`
- Definitions/theorems used:
  - `alternatingResolvent`
  - `AxisVanishesBelow`
  - `axisLinearOperator_increases_order`
  - `naturalOperator`
  - `naturalResolvent`
  - `naturalResolvent_equation`

The official resolvent is the alternating series `sum_k (-Q)^k`, where `Q = (1/2) J_2 M_chi`. The pinned filtration theorem states that one application of `Q` raises the vanishing radial order by one. Consequently, for radial row `n`, every term `Q^k A` with `k > n` vanishes at that row. Therefore the infinite series has the exact coordinate identity

`R(A)[n,m](eta) = sum_{k=0}^n (-1)^k (Q^k A)[n,m](eta)`.

The implementation uses this exact finite identity. It does not choose a numerical series tolerance or an arbitrary iteration cap.

## Actual-data anchoring

`actual_schedule_natural_resolvent(reference)` accepts only `ActualScheduleReferenceAxisState`. It reuses the already landed actual-schedule `AxisCoefficientNaturalOperator`, including the theorem-selected coefficient epsilon and the actual `chi` jets recovered from the pinned reference identity. Callers cannot inject a replacement `chi`, epsilon, operator table, or truncation threshold.

## Independent regression boundary

The regression suite checks:

1. the exact radial-filtration term count `n+1` for row `n`;
2. `naturalResolvent(one)` against the independently landed closed-form actual angular reference jets through several radial rows and parameter derivatives;
3. the inverse equation `R(A) + Q(R(A)) = A` on the actual axial reference state;
4. exact row-zero identity and fail-closed epsilon mismatch behavior;
5. explicit truth flags preserving `paper_exact=False` and the missing global AxisSpace certificate.

## Remaining boundary

This increment does **not** establish the global all-index weighted `AxisSpace` supremum/norm certificate or the Lean continuous-linear-map norm object. It also does not yet assemble the pinned `naturalRemainder`, evaluate `naturalRemainder(x0)`, produce a genuine Picard `x1`/fixed point, derive the final average/pressure fields, connect `NaturalProfileAssembly`, or verify the final support/moment/matching/cone obligations. `paper_exact_velocity_available` therefore remains false.
