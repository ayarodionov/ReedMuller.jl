# ReedMuller.jl

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

Candidates worth adding:
Sidel'nikov-Pershakov / derivative decoding for RM(2,m), permutation
(automorphism-group) decoding, Reed decoding with soft votes,
successive-cancellation viewing RM as polar codes with a different
frozen set, minimum-weight-parity-check ML for tiny codes.

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
* E. Abbe, A. Shpilka, M. Ye, "Reed-Muller codes: theory and
  algorithms", IEEE Trans. Inf. Theory, 2021 (survey).
