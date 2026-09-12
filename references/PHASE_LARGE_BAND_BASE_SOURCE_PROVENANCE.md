# Large-band base-source binding provenance

Status: **formal-structure**. This increment binds an already-admitted large-band LocalBase theorem witness to a stable sign-independent base-source identity. It does not materialize or certify the manuscript Proposition 5.5 background.

## Why this layer exists

`slow_labels.py` already encodes an important Section 6→7 dependency correctly: `TangentialBaseJetProvider.tangential_jet` is queried only from `(chart,R,Z,T)`. The label sign `sigma` does not enter the base-field provider; it is introduced later when the frozen representative phase is assembled. Thus the two wave labels `(ell,a,-1)` and `(ell,a,+1)` for one slow box are supposed to share one frozen base source.

The underflow-safe large-band path added in `phase_large_band_local_base.py` can represent `(ell,a,sigma)` and theorem-certified LocalBase ratios far beyond the binary64 chart cap, but before this increment it had no stable provider/source identity. A future caller could therefore attach the two signs of the same asymptotic slow box to unrelated background revisions while still passing each LocalBase admission separately.

`src/openai_ns_reconstruction/phase_large_band_base_source.py` closes only that identity gap.

## Interface and checks

`LargeBandBaseSourceWitness` records:

- the sign-free slow-box key `(ell,a)`;
- the exact LocalBase constant `M`;
- a stable `(source_id, source_revision)` key;
- theorem evidence kind and nonempty provenance;
- separate theorem assertions that the source obeys the normalized base-provider contract, is independent of `sigma`, and is exactly the source to which the admitted LocalBase bounds apply.

Only `analytic-theorem` and `formal-theorem` evidence tags are admitted. Sampled, fitted and numeric-scan evidence is rejected.

`LargeBandBaseSourceBinding` then requires exact agreement with an existing `LargeBandLocalBaseAdmission` in box key, band and `M`, rechecks that the LocalBase admission remains intact, and reconstructs the only two allowed signed labels from the common sign-free box key. The source revision is not inferred from field values and cannot be silently changed by the sign choice.

## Paper / repository mapping

Pinned paper/Lean context is unchanged from the preceding large-band LocalBase layer: `Finite Time Blowup for Navier-Stokes`, Sections 6–7.1 and the pinned `PhaseEstimates.LocalBaseBounds` contract. The executable repository-side sign-free provider convention is in `src/openai_ns_reconstruction/slow_labels.py`, where `TangentialBaseJetProvider.tangential_jet` accepts `chart,R,Z,T` and no `sigma` argument; `sigma` enters only when `RepresentativePhaseData` is created after freezing the jet.

This module does **not** assert that the current Stage 1/2 background reconstruction is already the Proposition 5.5 source. Its source witness is an admission/provenance contract for the future real source theorem.

## Independent regression

`tests/test_phase_large_band_base_source.py` uses formal fixtures only. It verifies that:

- one source binding yields exactly the `sigma=-1,+1` label pair with one common box/source key;
- the result is identical whether the pre-existing LocalBase witness happened to carry the plus or minus label;
- mismatched box, band or `M` is rejected;
- sampled/fitted/numeric-scan evidence, blank source/revision/provenance and any missing source-identity theorem fact are rejected;
- non-paper signs are rejected;
- no truth flag is upgraded by successful admission.

No test value is a sample from the manuscript flow.

## Remaining boundary

Still required before the Section 7 field estimates become paper-exact:

- a genuine Proposition 5.5/background provider with a pinned source revision;
- a machine-linked theorem showing that provider satisfies the admitted LocalBase ratios on every active enlarged slow box;
- the actual asymptotic Section 6 partition, representatives and support geometry;
- field-level vector/normal hypotheses for the frame and damping estimates;
- stress amplitudes, curl waves, mean corrections and Section 9 corrections built from those actual fields.

Therefore `actual_base_fields_verified=false`, `uniform_eq_7_9_to_7_11_verified=false`, and `paper_exact_velocity_available=false` remain mandatory.
