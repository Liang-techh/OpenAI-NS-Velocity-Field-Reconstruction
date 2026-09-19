# Repository replacement record

The user explicitly requested replacing the entire current contents of `Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction` with the strongest reproducible result, including visualizations, and making the presentation English.

The old main was observed at `065b67e9a3e8e22d697b49ccc67163db60a79783`. Before preparing the replacement, a preservation branch was created:

```text
archive/pre-st054-release-20260919
```

The replacement is a normal descendant commit, not an orphan history or force-push. The release branch is `release/st054-english-rebuild`. Old experimental files are removed from the current release tree rather than mixed with the new API, but remain recoverable from Git history and the archive branch. Historical issues, pull requests and research branches are not deleted or marked scientifically complete. The source repository with trailing `1` is not replaced.

## Selection

ST054-Q2 is the primary packaged checkpoint because it has the lower recorded volume L2 across the two controlled samples. ST054-M3 is retained because its sampled maximum is lower. Neither dominates every metric. ST055 is excluded until its complete field and evidence can be recovered and reproduced. Publication is not an optimization round.

## Included scope

Frozen velocity/pressure/force arrays, a small installable Python runtime, original independent validation, four replay reports, native MATLAB visualization/data, corrected dense-core tools, software tests and English documentation. The old repository's many partial paper-construction layers and work-in-progress notes are historical, not dependencies of the release. They are not silently rebranded as numerical success.

## Description and final publication receipt

`ABOUT.txt` contains the requested English repository description. Updating GitHub's repository About field is a separate repository-administration operation from committing README text. The final publication receipt must distinguish those two outcomes, along with the exact merged commit and actual completed CI jobs. No successful metadata update is inferred merely from this file existing.
