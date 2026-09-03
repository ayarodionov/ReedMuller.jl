---
name: run-tests
description: Run the ReedMuller.jl test suite and report the result. Use whenever the user asks to run tests, check that the package still passes, or verify that a change (new algorithm, refactor, doc edit) didn't break anything.
---

Run the full test suite from the package root:

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

Then:

- If every test passes, say so in one line (e.g. "All N tests pass") —
  no further explanation needed.
- If anything fails, show the specific failing `@testset` name(s) and
  the actual error/assertion output verbatim. Do not attempt to fix
  the failure unless the user asks you to; report it and stop.
- Do not add, remove, or edit any test as part of just running this
  skill.

If the user also wants the documentation build checked (this catches
a different class of break — e.g. a missing docstring for a newly
exported binding, which `checkdocs = :exports` in `docs/make.jl`
enforces), also run:

```
julia --project=docs docs/make.jl
```

from the package root, and report whether it completed without an
`ERROR:` line (a `Warning: Documenter could not auto-detect the
building environment. Skipping deployment.` line is expected locally
and is not a failure).
