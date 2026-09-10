# Repository working rules

1. Read `docs/NEXT_TASKS.md` and the machine-readable status before changing scope.
2. Never identify the bundled Gaussian example as OpenAI's constructed field.
3. Use independently specified forces/solutions for closure tests. Self-derived
   residual subtraction is not independent evidence.
4. Preserve `grad(c) cross A`, q-cutoff derivatives, axis regularity and the
   direct swirl's axisymmetry conditions.
5. New mathematical constructors need actual paper choices, provenance and
   support/moment/convergence evidence, not only interfaces or passing tests.
6. Run `python -m pytest -q -W error` and `python -m openai_ns_reconstruction verify`.
   The `--require-paper-exact` gate must keep failing until the real complete
   construction is instantiated and its evidence is reviewed.
7. Distinguish local tests from remote CI and a local patch from a pushed commit.
   Do not overwrite concurrent work or use force-push to integrate this patch.
