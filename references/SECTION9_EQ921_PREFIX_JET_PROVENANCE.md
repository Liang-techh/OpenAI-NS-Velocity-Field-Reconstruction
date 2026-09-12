# Section 9 Eq. (9.21) finite-prefix analytic-jet provenance

Status: **formal-structure / analytic adapter only**. This increment does not
materialize the manuscript's actual Section 7/8 correction fields and does not
make the paper-exact velocity available.

## Scope

`section9_eq921_prefix_jet.py` bridges the already-landed point-value finite
prefix into a derivative-aware representation. For each admitted positive
stage it accepts a complete provider-supplied jet of the weighted Eq. (9.21)
summand

`chi(a_j q) (A_j, B_j, p_j)`

with derivatives indexed by physical `(t,x,y,z)` multi-indices through one
caller-selected total order. The adapter performs no finite differences, fits,
or sampled derivative reconstruction.

The bridge fails closed unless:

1. stages cover exactly the same contiguous prefix as the existing
   `Section9FinitePrefixEvaluationCertificate`;
2. each stage scale and exact `q` agree with that evaluation and the existing
   correction-stage admission;
3. the `(A,B,p)` evidence kinds and provenance strings are identical to the
   corresponding admission certificate;
4. every stage uses one common provider id/revision/provenance and one common
   derivative order;
5. all `(dt,dx,dy,dz)` multi-indices through that order are present; and
6. each weighted stage jet's zero-order `(A,B,p)` value agrees exactly with the
   independently landed finite-prefix contribution before any derivatives are
   summed.

The resulting prefix jet is then a literal termwise sum of the supplied
analytic derivative tables.

## Why this is not the residual yet

The higher derivative entries remain provider inputs. This module does not
machine-derive them from the actual Section 7/8 correction construction, does
not prove that the provider differentiated the manuscript's fixed cutoff, and
does not construct the locally finite infinite Eq. (9.21) field through the
endpoint. Consequently it is not valid to feed a synthetic fixture through
this adapter and label the output the paper residual.

The following truth boundaries remain false:

- `analytic_derivatives_machine_derived_from_actual_corrections`
- `paper_fixed_cutoff_derivatives_machine_verified`
- `actual_correction_field_values_verified`
- `actual_section9_sequence_verified`
- `eq_9_21_infinite_sum_constructed`
- `endpoint_covered`
- `residual_artifact_ready`
- `paper_exact_velocity_available`

## Next paper-exact dependency

A real Section 7/8 correction-field provider must export the actual weighted
summand values and analytic spacetime derivatives, with the same stage/source
identity already admitted by Section 9. Once that provider populates this
adapter, a later step can assemble the required local velocity/pressure jets
and apply the Navier--Stokes operator to produce a genuine residual artifact.
Only after that can the existing Eq. (9.18) envelope/endpoint-majorant path be
fed by machine-derived residual data rather than theorem sidecars or fixtures.
