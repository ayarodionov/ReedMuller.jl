---
name: add-algorithm
description: Add a new Reed-Muller encoding or decoding algorithm to this package, following its established conventions (common interface, tests, docs, README). Use when the user asks to implement/add a new encoder, decoder, or generic decoder wrapper — from a paper, a PDF, or a plain description.
argument-hint: [algorithm name or description, and a paper/PDF reference if you have one]
---

Add the algorithm described in `$ARGUMENTS` to ReedMuller.jl, following
every convention below. Treat this as a checklist, not a suggestion —
every past algorithm in this package (Reed, FHT, Dumer, Dumer-Shabunov,
Sidel'nikov-Pershakov, RPA, BP, GLP, Kamenev's graph search, and the
generic wrappers) was added this exact way, and the test/docs/CI setup
assumes it.

## 1. Understand it before writing code

If given a PDF, read it in full with the Read tool (it handles PDFs
directly — extract by page range for figures/algorithm boxes rather
than guessing). Work out precisely: the message-coordinate convention
it expects or produces (`:monomial` — coefficients of Boolean
monomials, degree ascending — or `:plotkin` — recursive (u, u+v)
layout, u-part first), whether it's soft-input (LLRs) or
hard-decision, its parameters and their defaults, and any worked
example in the paper you can reproduce as a correctness test.

**Verify the citation independently before writing it anywhere.** A
web search's AI-generated summary can be subtly wrong even when it
sounds authoritative (this happened in this repo: a search summary
gave the wrong page range for the Sidel'nikov-Pershakov paper, twice,
from two different searches — it was only caught by reading the
actual reference list of a real paper that cites it). If you can
fetch or already have the actual paper, read its own header/DOI. If
not, cross-check the citation against another paper's real
bibliography (WebFetch it and read the references section), not
just a second search summary.

## 2. Place the code

- A standalone decoding algorithm → new file `src/decoders/<name>.jl`.
- A generic wrapper that composes with *any* inner decoder (it takes
  `code` and an `inner::AbstractDecoder`/`AbstractEncoder` argument
  rather than implementing decoding itself) → add to
  `src/decoders/wrappers.jl` if small, or its own file if substantial
  (see `src/decoders/aed.jl` or `src/decoders/glp.jl` for examples of
  the latter).
- A new encoder → `src/encoders.jl`.

## 3. Implement the common interface

Every decoder is a `struct ... <: AbstractDecoder` (encoders:
`AbstractEncoder`), defined in `src/code.jl`, with:

- An inner or keyword constructor that validates every parameter and
  throws `ArgumentError` with a clear message for invalid values
  (see any existing decoder's constructor for the pattern:
  `x in (:a, :b) || throw(ArgumentError("..."))`).
- `basis(::T) = :monomial` or `:plotkin` — the message convention this
  decoder's `decode` expects.
- `decode(dec::T, code::RMCode, llr::AbstractVector{<:Real}) -> BitVector`
  (or `encode(enc::T, code::RMCode, message::AbstractVector) -> BitVector`).
  First line inside: validate the input length and throw
  `DimensionMismatch` on mismatch — copy the exact pattern from any
  existing decoder (`n = blocklength(code); length(llr) == n || throw(DimensionMismatch(...))`).
- A full docstring above the struct: the constructor signature, one
  or two paragraphs describing the algorithm (what it does, not just
  what it's called), the verified paper citation, its complexity, and
  a final `Message convention: `:monomial`.` (or `:plotkin`) line.
  Look at `src/decoders/sidelnikov_pershakov.jl` or
  `src/decoders/graph_search.jl` for the level of detail expected.

LLR sign convention (already fixed by the package, don't reinvent
it): positive LLR = bit 0 more likely. `hard_llr(bits)` converts plain
bits to pseudo-LLRs for decoders that only need hard input.

## 4. Wire it into the module

In `src/ReedMuller.jl`:
- Add `include("decoders/<file>.jl")` (or the encoder file) in a
  sensible position relative to what it depends on.
- Add every new exported type/function to the `export` list.

## 5. Tests — `test/runtests.jl`

- Add `(code, encoder, decoder)` tuples to the shared noiseless
  roundtrip `cases` list, covering at least two code sizes and, where
  meaningful, an edge case (`r = 0` or `r = m`).
- Add a dedicated `@testset "<algorithm> decoding"` with
  behavior-specific checks, not just "it runs":
  - If the paper has a worked example, reproduce it exactly (see the
    `"graph search decoding"` testset — it decodes the paper's own
    example LLR vector and checks the exact resulting codeword).
  - Compare against `MLDecoder` on a code small enough for brute
    force (k ≤ 24), or against an existing near-ML decoder on a noisy
    batch, to confirm it actually performs as claimed (not just that
    it doesn't crash) — see the `RPA decoding` or `GLP decoding`
    testsets for the pattern (run N trials at a fixed Eb/N0 via
    `simulate` or a manual loop, compare error counts).
  - Add `@test_throws ArgumentError ...` for every invalid
    constructor argument and `@test_throws DimensionMismatch
    decode(...)` for a wrong-length LLR vector.
- Run the `run-tests` skill (or
  `julia --project=. -e 'using Pkg; Pkg.test()'` directly) and confirm
  every test passes before moving on.

## 6. Documentation

- `ALGORITHMS.md`: add a `### <Algorithm name>` section under the
  right heading (Encoders / Decoders / Generic decoder wrappers),
  matching the existing format exactly: a `**File:**` line linking to
  the source file, a `**Reference:**` line with the verified
  citation, then the same prose description as the docstring.
- `docs/src/api.md`: add every new exported name to the appropriate
  `` ```@docs ``` `` block. This is not optional — `docs/make.jl` sets
  `checkdocs = :exports`, so a missing entry fails the Documentation
  CI workflow, not just a style nit.
- `README.md`: add a row to the "Implemented algorithms" table
  (`| Algorithm | Type | Input | Basis | Scope | Complexity |`) and, if
  it fills a gap, remove it from the "Candidates not yet implemented"
  list at the end of `ALGORITHMS.md` and the README if it's listed
  there.
- Verify the docs build locally before committing:
  `julia --project=docs docs/make.jl` — must complete with no `ERROR:`
  line (the `Documenter could not auto-detect the building
  environment. Skipping deployment.` warning is expected and fine
  locally).

## 7. Benchmarks (only if the user asks, or it's clearly comparison-worthy)

- `benchmarks/compare.jl` for the small default sweep — add a
  `(name, encoder, decoder)` tuple to its `pipelines` list.
- `benchmarks/RESULTS.md` and `benchmarks/plots/generate_charts.jl` are
  a much larger undertaking (100,000 trials/point across RM(2,8-10),
  can take 30+ minutes and has previously hit machine memory limits —
  check `sysctl kern.memorystatus_vm_pressure_level` before launching
  and use a `kill -0 $pid` watchdog for any long background run since
  a killed process leaves no error text). Don't add to these unless
  explicitly asked; if asked, follow the pattern of the existing
  `benchmarks/compare_large_*.jl` scripts.

## 8. Before calling it done

1. All tests pass (step 5).
2. Docs build cleanly (step 6).
3. `git status` shows only the files you intended to touch.
4. Commit with a message describing what was added and how it was
   verified (see `git log` in this repo for the established style —
   summary line, blank line, a short paragraph on what was
   implemented and what tests/checks confirm it works). Only push if
   that's the established workflow in this conversation/repo, or the
   user asks.
