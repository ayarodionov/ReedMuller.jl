# ReedMuller.jl

[![CI](https://github.com/ayarodionov/ReedMuller.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ayarodionov/ReedMuller.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://ayarodionov.github.io/ReedMuller.jl/)

A common place — and a common language — for Reed-Muller encoding and
decoding algorithms in Julia. The goal is to collect as many
implementations as possible behind one interface so they can be
compared head-to-head under identical conditions.

## The common interface

```julia
using ReedMuller

code = RMCode(2, 5)                 # RM(r, m): n = 2^m, k = Σ C(m,i), d = 2^(m-r)

enc = MatrixEncoder(code)           # any AbstractEncoder
msg = rand(Bool, dimension(code))
cw  = encode(enc, code, msg)        # BitVector of length 2^m

ch  = BIAWGN(0.8)                   # or BSC(p)
llr = transmit(Random.default_rng(), ch, cw)

dec = ReedDecoder()                 # any AbstractDecoder
est = decode(dec, code, llr)        # recovered message bits
```

All decoders take a vector of **log-likelihood ratios** (positive =
bit 0 more likely). Hard-decision decoders threshold internally; use
`hard_llr(bits)` to feed a plain bit vector to any decoder.

### Message bases

Two message-coordinate conventions coexist, reported by `basis(x)`:

* `:monomial` — message bits are coefficients of Boolean monomials,
  degree ascending (Reed, FHT decoders).
* `:plotkin` — recursive (u, u+v) layout, u-part first (Dumer decoder).

Both bases generate the same code; only the message-to-codeword map
differs. `simulate` refuses to pair an encoder and decoder with
mismatched bases.

## Implemented algorithms

| Algorithm | Type | Input | Basis | Scope | Complexity |
|---|---|---|---|---|---|
| `MatrixEncoder` | encoder | — | either | any RM(r,m) | O(nk) |
| `PlotkinEncoder` | encoder | — | `:plotkin` | any RM(r,m) | O(n log n) |
| `ReedDecoder` | decoder | hard | `:monomial` | any RM(r,m) | O(nk) |
| `FHTDecoder` | decoder | soft, ML | `:monomial` | RM(1,m) only | O(n log n) |
| `DumerDecoder` | decoder | soft | `:plotkin` | any RM(r,m) | O(n log n) |
| `DumerShabunovDecoder` | decoder | soft, list | `:plotkin` | any RM(r,m) | O(L·n log n) |
| `SidelnikovPershakovDecoder` | decoder | soft | `:monomial` | RM(2,m) only | O(n² log n) |
| `RPADecoder` | decoder | soft, near-ML | `:monomial` | r ≥ 1 (best r ≤ 3) | O(it·n² log n)/level |
| `BPDecoder` | decoder | soft, iterative | `:monomial` | any RM(r,m) | O(it·Σ row wt) |
| `MLDecoder` | decoder | soft, exact ML | either | k ≤ 24 | O(2^k·n) |
| `GLPDecoder` | decoder | soft, list+perms | `:plotkin` | any RM(r,m) | O(L·n log n) |

Options: `DumerDecoder`/`DumerShabunovDecoder` take `leaves = :fht` to
terminate the recursion at first-order nodes with (list-)FHT decoding
— the full Dumer-Shabunov construction, noticeably stronger;
`SidelnikovPershakovDecoder` takes `voting = :majority` for Sakkour's
simplified plain-majority variant.

`GLPDecoder(code, L; perms)` is the global-list-with-permutations
decoder after the `dtrm_glp` codec of
[ecclab](https://github.com/kshabunov/ecclab): the Dumer-Shabunov
list is seeded with one path per code-preserving variable permutation
of the received word, and all paths then compete for the same `L`
slots. `perms` takes ecclab's ensembles — `:identity` (`P1`,
equivalent to plain Dumer-Shabunov), `:cyclic` (`Pcyclic`), `:pairs`
(`PmCr`, every variable pair moved to the decided-first positions) —
or an explicit list of variable permutations. Ecclab's decoders map
to this package as: `dtrm0` ≈ `DumerDecoder()`, `dtrm1` ≈
`DumerDecoder(leaves = :fht)`, `dtrm_glp` ≈ `GLPDecoder`, `rm1_ml` ≈
`FHTDecoder`. Note: the ensembles pay off when `L` is well above the
ensemble size (ecclab runs `PmCr` with L = 256 at m = 8-9); with a
tight budget the shared list dilutes and plain Dumer-Shabunov can win.

### Generic wrappers

These compose with *any* decoder through the common interface:

| Wrapper | Idea | Cost |
|---|---|---|
| `AutomorphismEnsembleDecoder(code, inner; size)` | decode `size` permuted received words (automorphism group GA(m)), keep best by correlation | size × inner |
| `ChaseDecoder(code, inner; t)` | try all 2^t sign flips of the t least reliable positions | 2^t × inner |
| `GMDDecoder(code, inner)` | erase 0, 2, 4, … least reliable positions, keep best | (d+1)/2 × inner |

## Comparing algorithms

```julia
res = simulate(PlotkinEncoder(), DumerDecoder(), code, BSC(0.02); trials = 100_000)
# res.ber, res.wer
```

`benchmarks/compare.jl` sweeps Eb/N0 over the BI-AWGN channel and
prints a WER table for every pipeline:

```
julia --project=. benchmarks/compare.jl
```

## Adding a new algorithm

1. Create `src/decoders/<name>.jl` (or add an encoder in
   `src/encoders.jl`).
2. Define a struct subtyping `AbstractDecoder` / `AbstractEncoder`
   with any algorithm parameters as fields.
3. Implement `decode(dec, code, llr)::BitVector` (or
   `encode(enc, code, msg)::BitVector`) and `basis(x)`.
4. `include` the file in `src/ReedMuller.jl`, export the type, and add
   a noiseless-roundtrip case in `test/runtests.jl` plus a pipeline in
   `benchmarks/compare.jl`.

Candidates worth adding: Kabatiansky-Tavernier (and Fourquet-Tavernier)
deterministic list decoding of RM(2,m), sparse/multi-decoder RPA and
the hardware-oriented IPA variant, Reed decoding with soft votes,
Sidel'nikov-Pershakov with candidate lists per derivative.

## Running the tests

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

## References

* I. S. Reed, "A class of multiple-error-correcting codes and the
  decoding scheme", 1954.
* R. R. Green, "A serial orthogonal decoder", JPL Space Programs
  Summary, 1966 (FHT decoding of first-order codes).
* I. Dumer, "Recursive decoding and its performance for low-rate
  Reed-Muller codes", IEEE Trans. Inf. Theory, 2004.
* I. Dumer, K. Shabunov, "Soft-decision decoding of Reed-Muller
  codes: recursive lists", IEEE Trans. Inf. Theory, 2006.
* K. Shabunov, ecclab — recursive decoding simulation programs
  (reference C implementations of dtrm0/dtrm1/dtrm_glp),
  https://github.com/kshabunov/ecclab.
* V. M. Sidel'nikov, A. S. Pershakov, "Decoding of Reed-Muller codes
  with a large number of errors", Probl. Inf. Transm., 1992.
* B. Sakkour, "Decoding of second order Reed-Muller codes with a
  large number of errors", IEEE ITW, 2005.
* M. Ye, E. Abbe, "Recursive projection-aggregation decoding of
  Reed-Muller codes", IEEE Trans. Inf. Theory, 2020.
* M. Geiselhart, A. Elkelesh, M. Ebada, S. Cammerer, S. ten Brink,
  "Automorphism ensemble decoding of Reed-Muller codes", IEEE Trans.
  Commun., 2021.
* D. Chase, "A class of algorithms for decoding block codes with
  channel measurement information", IEEE Trans. Inf. Theory, 1972.
* G. D. Forney, "Generalized minimum distance decoding", IEEE Trans.
  Inf. Theory, 1966.
* E. Abbe, A. Shpilka, M. Ye, "Reed-Muller codes: theory and
  algorithms", IEEE Trans. Inf. Theory, 2021 (survey).
