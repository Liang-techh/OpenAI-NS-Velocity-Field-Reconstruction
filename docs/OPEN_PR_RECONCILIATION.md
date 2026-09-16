# Open PR reconciliation at checkpoint `77ed17c`

Task: NS020  
Owner: Agent 4  
Checkpoint branch: `codex/stage1-complete-slow2`  
Checkpoint SHA: `77ed17c1e81b7abce99fa4c7ffc5b3d96ba389d0`

This is a read-only reconciliation of the legacy/open PR graph against the published checkpoint. It does **not** accept, merge, close, or promote any mathematical status. `paper_exact_velocity_available=false` and `full_reconstruction=false` remain unchanged.

## Status meanings used here

- **implementation present**: the PR contains the advertised narrow code/documentation payload.
- **coordinator acceptance**: independent acceptance under the new #368 protocol. None of the legacy PRs below has that acceptance merely because CI passed.
- **merge status**: whether the PR is merged. Every priority PR below is still open/unmerged.
- A successful old CI run verifies that exact old head only. It does not certify the checkpoint branch or a future replay.

## Priority open PRs

| PR | Scope | Base -> head | Exact-head CI observed | Checkpoint comparison | Local delivery / acceptance / merge |
| --- | --- | --- | --- | --- | --- |
| #366 | Stage-1 genuine `slow1(x1)`, stacked after the first-Picard constituent chain | `agent6-first-picard-lin2-x1` `f070402da90a454d12f944b62dcb10a4b919b345` -> `agent6-first-picard-slow1-x1` `a30a8ad8aba31724865d662411badb5c1a3e19b3` | run `34862205852`: **failure**. Both full Python/NumPy jobs failed in `Test full suite and dependencies`; slice jobs observed in the run include successes. | Diverged from checkpoint: 8 ahead / 18 behind, merge base `13f9c53b...`. Comparison carries six first-Picard constituent source/test/provenance families not present on the checkpoint branch. | implementation present; **not coordinator-accepted**; open/unmerged. Do not integrate while full exact-head CI is red. |
| #365 | Stage-2 hierarchy-owned third-eta PositiveAxis forcing replay | `main` `065b67e9a3e8e22d697b49ccc67163db60a79783` -> `agent2/replay-forcing-third-eta-main-065b` `7d47e3032467ff5302d055dc4ac478ce706e7e15` | run `34861426966`: **success** | Diverged: 4 ahead / 10 behind checkpoint; merge base `065b67e...`. Four-file delta: new forcing-third-eta source/provenance plus extension of the lower-source-third-eta regression. | narrow implementation present and green on its head; acceptance pending; open/unmerged. Must be replayed onto the checkpoint, not merged as a stale-main snapshot. |
| #364 | Stage-1 exact-average axis callback replay | `main` `065b67e9a3e8e22d697b49ccc67163db60a79783` -> `agent1/exact-average-axis-callback-065b67e` `f0ff5b9717bf5a470f31e6e0cae78af55cf4b6cf` | run `34859620484`: **success** | Diverged: 1 ahead / 10 behind; merge base `065b67e...`. It modifies existing checkpoint file `src/openai_ns_reconstruction/profiles.py` and adds one focused test plus provenance. | implementation present/green; acceptance pending; open/unmerged. Relevant later to NS016; replay manually after the genuine fixed-point/profile path exists. |
| #363 | Section-10 all-order global small-scale residual majorant, stacked on #360 | `agent4/exterior-residual-zero-germ-d187be8` `f2f7822c4df9c78ec645ed88908b37ad76fefaf4` -> `agent4/global-endpoint-residual-majorant-f2f7822` `c52e0193017603dbe4bae159f4b2505439b551cf` | run `34858962367`: **success** | Diverged: 12 ahead / 11 behind; merge base `1bc98e6...`. Relative to checkpoint the stack adds the PaperLocalization spine, exterior zero-germ, and global endpoint-majorant modules/tests/provenance. | theorem-shaped implementation present/green; acceptance pending; open/unmerged. It is not runtime residual/force closure. |
| #362 | Agent-7 finite all-jet physical target boxes, stacked on the cofinal finite-tail ladder | `agent7/cofinal-physical-tail-ladder` `68f41644b6b9dbb992f9ef36bfae2894e22e9a7f` -> `agent7/finite-physical-tail-target-box` `1d66573cb73f96c83ccbbeca8219abb17cff1d8d` | run `34856682255`: **cancelled** | Diverged: 60 ahead / 21 behind; merge base `51e7b91...`. The comparison contains the long finite target/cutoff/majorant stack; it is not an infinite provider. | implementation present, but no successful exact-head run from the observed run; acceptance pending; open/unmerged. Do not merge the long stack wholesale. |
| #361 | Prepared/LocalBase formal binding replay | `main` `065b67e9a3e8e22d697b49ccc67163db60a79783` -> `agent5/prepared-localbase-replay-065b67e` `8ad6ed38352072ad8849b0f42e2fba82dd9277f3` | run `34856228986`: **success** | Diverged: 4 ahead / 10 behind; merge base `065b67e...`. Four unique Prepared LocalBase binding/provenance/test files are absent from checkpoint. | implementation present/green; acceptance pending; open/unmerged. It remains a formal binding, not materialized LocalBase values/jets. |
| #360 | Section-10 exterior residual zero-germ, stacked on the PaperLocalization spine (#355) | `agent4/paper-localization-spine-1bc98e6` `d187be8783d5a1aba894d03b9ea14f6276e14994` -> `agent4/exterior-residual-zero-germ-d187be8` `f2f7822c4df9c78ec645ed88908b37ad76fefaf4` | run `34852698191`: **success** | Diverged: 6 ahead / 11 behind; merge base `1bc98e6...`. Relative to checkpoint the stack adds the PaperLocalization spine and zero-germ layers. | implementation present/green; acceptance pending; open/unmerged. It is an ancestor of #363; do not integrate #360 and #363 independently. |

## Replay / supersession relationships

- #365 is the current-main replay replacement for closed/unmerged #322. Use #365, not #322, if this narrow forcing seam is later consumed.
- #364 is the current-main replay replacement for closed/unmerged #356. Use #364, not #356.
- #361 is the current-main replay replacement for closed/unmerged #359 (itself a replay of earlier Prepared LocalBase work). Use #361, not #359.
- #363 extends the #355 -> #360 Section-10 stack. If that formal chain is later replayed, use the top required stack content once; do not merge #360 and then separately replay the same ancestor through #363.
- #366 is a stacked Stage-1 chain whose direct base is #357; the chain also includes the earlier source/pressure/quad1/lin1/lin2 increments. Its top head diverges substantially from the checkpoint and currently has failing full CI.
- #362 is the top of a long finite-target/cutoff/majorant stack. It remains finite/conditional and must not be treated as NS025/NS026 all-order closure.

## New #368 task-board PRs are a separate integration lane

At this reconciliation point the new task-board work targets the checkpoint directly:

- #369 / NS001: base `77ed17c`, head `376a6e551a9aa373c0f0ad5fe55b7c34c54d0852`, open.
- #370 / NS003: base `77ed17c`, head `7309bef0f569d7a2b3149afe40dfc24d8784808e`, open; its final-head full CI was still running when inspected.
- #371 / NS006: base `77ed17c`, head `5f63f32b934d8399432f858f7dd62cb7b7c03ff6`, open; focused local report says 15 passed, while full-repository checks were not locally run.

These checkpoint-targeted PRs should not be overwritten by merging old main-based or stacked branches.

## Minimal integration order

1. **Keep `codex/stage1-complete-slow2` as the sole integration baseline.** Do not merge legacy main-based/stacked PRs directly if doing so would discard the checkpoint's 10-21 newer commits.
2. **Finish and coordinator-review the active checkpoint-targeted NS001/NS003/NS006 work first.** These are on the live dependency chain for strict Stage-1 arithmetic.
3. **Do not integrate #366 now.** Its exact-head full suites failed and its long stack is far behind the checkpoint. When NS013/NS015 reaches the first-Picard/fixed-point seam, replay only the still-needed constituent implementations onto the then-current checkpoint, resolve against the strict interval data, and require fresh full CI.
4. **Hold #364 for NS016.** It is a small identity fix with green old-head CI, but the meaningful consumer is the eventual genuine profile assembly. Replay its `profiles.py` change then, rather than merging a stale-main snapshot now.
5. **Hold #365 for NS022.** It is a narrow green four-file Stage-2 replay and is the preferred replacement for #322, but NS022 depends on Stage-1 closure (NS017). Replay it onto the live checkpoint when that dependency is actually satisfied.
6. **Treat #362 as reference material for NS025/NS026, not as all-order completion.** Its current observed run was cancelled and the stack is 60 commits ahead / 21 behind the checkpoint. Reuse the finite certificates only after an actual hierarchy-owned total provider exists.
7. **Hold #361 for NS027.** Its formal Prepared/LocalBase bridge is clean and green, but real background values/jets first require NS026. Replay the four-file payload at that stage.
8. **Hold the #355 -> #360 -> #363 chain for NS034.** #363 is the top useful formal implication and its old head is green, but runtime force/endpoint closure depends on the actual NS033 field. Replay the required chain once onto the then-current checkpoint; do not merge ancestor and descendant redundantly.
9. **Run NS021 only on an explicitly selected integration commit.** The checkpoint's own run `35045650971` was still in progress when inspected; no old PR run may be reported as validation of a different integrated tree.

## Conflict / safety notes

- A green PR head is evidence only for that exact head; coordinator acceptance remains pending until replay/integration evidence is checked.
- A cancelled run (#362) is neither success nor code failure; it simply provides no successful exact-head acceptance evidence.
- A failing run (#366) blocks integration until the failing full suites are understood and a fresh head is green.
- The checkpoint comparison shows every priority legacy head has diverged from `77ed17c`; none should be merged as though it were already based on the checkpoint.
- Formal/theorem-shaped Section 5-10 adapters must retain their current truth boundaries. Replaying them does not by itself make actual fields, infinite sums, force artifacts, or paper-exact reconstruction available.
