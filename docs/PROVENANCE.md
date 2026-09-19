# Source and artifact provenance

The public research source is `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1`. This publication target is the separate repository without the trailing `1`.

| Component | Pinned source |
|---|---|
| ST054 study, recipes and numerical operator | `c77492a48e9c0f13d4d51244987c57c28519409b`, research PR #630 |
| Native MATLAB evaluator/viewer and tested MAT export | `5d49ccf07d097a4bffeb35f6398936493b062d00`, merged research PR #689 |
| Previous target repository main | `065b67e9a3e8e22d697b49ccc67163db60a79783` |

The complete path-level mapping and SHA-256 digests are in `tools/source_manifest.json`. `tools/build_release.py` verifies upstream bytes, applies only namespaced imports to the Python runtime and derives release JSON files from the frozen MAT coefficient arrays. It does not execute an optimizer. The resulting `evidence/import_receipt.json` records both source and release hashes after setup; it is already present in the initialized complete ZIP.

The canonical MAT SHA-256 is

```text
59d43fdc40cd4df1675c2191f85810b9351802b2f4b1e0ac25a3bccfd9914c62
```

Original research candidate JSON identities, retained as provenance:

```text
ST054-Q2 d772949621e0f7703a9d6b36f28828532ba3e5ec678dec14db5ce14eb8b851d5
ST054-M3 2b2b571986297b52bc884a2ee4808a7ade1ade86f38f938604f2a64118882df3
```

Release JSON metadata and small ancestor-replay floating-point differences mean these are not the hashes of the new release JSON files. After one-time setup, the generated release files have their own hashes under `src/ns_reconstruction/data/manifest.json`; the initialized ZIP already contains them. Their entire coefficient arrays equal those stored in the pinned MAT file. Velocity, pressure and residual are compared with the export's reference arrays at `1e-8` absolute tolerance, and the four original scientific rejections are separately replayed.

Original MATLAB receipts are explicitly stored under `evidence/upstream_matlab/`, not relabeled as current release tests. New CI artifacts and any checked-in release receipts name their own run, commit and scope. The corrected core tools are new release code, not silently attributed to the old upstream tests.

The earlier full-paper/exact-reconstruction work remains in the target's Git history and archive branch. The release does not import unverified theorem completion flags or claim that formal adapters constitute an implemented original velocity field. The ongoing source experiment archive remains separate.

No third-party paper PDF, source numerical coefficient identity, external research solver execution, Lean build or new license grant is invented by this migration. Existing author attribution is retained in the source history. This repository's name is historical; the project is unaffiliated with OpenAI.
