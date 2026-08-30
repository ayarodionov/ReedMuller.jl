# Algorithms

Every encoder and decoder in this package implements the common
interface described in the [README](README.md): `encode`/`decode`
over a fixed `RMCode`, LLR input, and a declared message `basis`
(`:monomial` or `:plotkin`). This file gives a fuller description of
each algorithm than fits in the README's summary table, and points to
its implementation and originating paper.

## Encoders

### Matrix encoder
**File:** [`src/encoders.jl`](src/encoders.jl) — `MatrixEncoder`

Reference encoder: multiplies the message by a precomputed generator
matrix over F₂. Works with either message basis (`:monomial` or
`:plotkin`, selected at construction).

### Plotkin encoder
**File:** [`src/encoders.jl`](src/encoders.jl) — `PlotkinEncoder`

Fast recursive (u, u+v) encoder: O(n log n) XORs, no stored matrix.
Message convention is `:plotkin` (u-part bits first, then v-part),
identical to `MatrixEncoder(code; basis=:plotkin)`.

## Decoders

### Reed's majority-logic decoder
**File:** [`src/decoders/reed.jl`](src/decoders/reed.jl) — `ReedDecoder`
**Reference:** I. S. Reed, "A class of multiple-error-correcting codes
and the decoding scheme", 1954.

Hard-decision decoding: LLRs are thresholded to bits, then monomial
coefficients are recovered degree by degree, from degree r down to 0.
The coefficient of each degree-d monomial x_S is the majority over
2^(m-d) parity checks, one per assignment of the variables outside S;
after a degree is finished its contribution is subtracted from the
received word. Guaranteed to correct up to ⌊(2^(m-r) - 1) / 2⌋ errors.
Message convention: `:monomial`.

### Fast Hadamard transform (FHT) decoder
**File:** [`src/decoders/fht.jl`](src/decoders/fht.jl) — `FHTDecoder`
**Reference:** R. R. Green, "A serial orthogonal decoder", JPL Space
Programs Summary, 1966 (the "Green machine").

Maximum-likelihood soft-decision decoder for first-order codes
RM(1, m). The transform of the LLR vector is maximized in absolute
value at the index whose bits are the linear coefficients; the sign
gives the constant term. Complexity O(n log n). Message convention:
`:monomial`, i.e. `[a₀, a₁, …, a_m]` with codeword
`c(z) = a₀ ⊕ ⊕ᵢ aᵢ zᵢ`.

### Dumer's recursive decoder
**File:** [`src/decoders/dumer.jl`](src/decoders/dumer.jl) — `DumerDecoder`
**Reference:** I. Dumer, "Recursive decoding and its performance for
low-rate Reed-Muller codes", IEEE Trans. Inf. Theory, 2004.

Splits RM(r, m) via the Plotkin (u, u+v) decomposition, first decoding
the v-part from combined half-LLRs, then the u-part from the
v-corrected sum. Leaf nodes: repetition codes (r = 0) decide by the
LLR sum, rate-1 codes (r = m) by symbol-wise hard decision.
`combine` chooses the check-node rule for the v-branch (`:minsum` or
exact `:exact` LLR combining); `leaves = :fht` stops the recursion at
first-order nodes and decodes them optimally with the FHT instead of
recursing to single bits. Complexity O(n log n). Message convention:
`:plotkin`.

### Dumer-Shabunov recursive list decoder
**File:** [`src/decoders/dumer_shabunov.jl`](src/decoders/dumer_shabunov.jl) — `DumerShabunovDecoder`
**Reference:** I. Dumer, K. Shabunov, "Soft-decision decoding of
Reed-Muller codes: recursive lists", IEEE Trans. Inf. Theory, 2006.

Same Plotkin recursion as `DumerDecoder`, but instead of committing to
one decision per node it maintains a list of up to `L` candidate
decoding paths. Information decisions are made at the leaves —
repetition nodes and single-position nodes, to which rate-1 subcodes
reduce recursively — where every surviving path branches into both
hypotheses. Paths are ranked by the standard successive-cancellation
path metric (accumulated magnitude of contradicted leaf LLRs); the
list is pruned back to the `L` best paths after each branching. With
`leaves = :fht` the recursion instead terminates at first-order nodes,
expanding every path into the full list of affine codewords scored by
exact correlation (the paper's list-FHT termination). `L = 1` with
`:bits` leaves reduces to plain `DumerDecoder`; growing `L` approaches
ML. Complexity O(L·n·log n). Message convention: `:plotkin`.

### Sidel'nikov-Pershakov derivative decoder
**File:** [`src/decoders/sidelnikov_pershakov.jl`](src/decoders/sidelnikov_pershakov.jl) — `SidelnikovPershakovDecoder`
**References:** V. M. Sidel'nikov, A. S. Pershakov, "Decoding of
Reed-Muller codes with a large number of errors", Probl. Inf. Transm.,
1992; B. Sakkour, "Decoding of second order Reed-Muller codes with a
large number of errors", IEEE ITW, 2005 (the `:majority` voting mode).

For second-order codes RM(2, m) only. A codeword is a quadratic
Boolean function f(x) = xᵀQx ⊕ ⟨l, x⟩ ⊕ a₀, so its derivative in any
direction b, `D_b f(x) = f(x) ⊕ f(x ⊕ b) = ⟨Bb, x⟩ ⊕ f(b) ⊕ f(0)` with
`B = Q + Qᵀ`, is affine:

1. For every nonzero direction b, form derivative LLRs by min-sum
   combining and ML-decode them in RM(1, m) via FHT, giving an
   estimate `w_b ≈ Bb` with a reliability (the FHT peak magnitude).
2. Recover each entry `B_ij` by a majority vote over all pairs
   `(b, b ⊕ e_i)`, whose estimates XOR to `B e_i` — `voting = :weighted`
   scales each vote by the reliability of its less reliable member,
   `:majority` counts every vote equally (Sakkour's simplification).
3. Peel the recovered quadratic part off the channel LLRs and
   FHT-decode the residual RM(1, m) for the affine part.

Corrects any pattern of up to ⌊(d-1)/2⌋ = 2^(m-3) - 1 errors and, on
random errors, most patterns well beyond half the minimum distance.
Complexity O(n² log n). Message convention: `:monomial`.

### Recursive Projection-Aggregation (RPA) decoder
**File:** [`src/decoders/rpa.jl`](src/decoders/rpa.jl) — `RPADecoder`
**Reference:** M. Ye, E. Abbe, "Recursive projection-aggregation
decoding of Reed-Muller codes", IEEE Trans. Inf. Theory, 2020.

Near-ML for low orders (r ≤ 3). For RM(r, m) with r ≥ 2, each
iteration: **(1) Projection** — for every nonzero direction z, project
the received LLRs onto the quotient by the subspace {0, z} (the LLR of
`y(x) ⊕ y(x ⊕ z)` per coset), which carries a codeword of
RM(r-1, m-1); **(2) Recursion** — decode each projection recursively,
reaching RM(1, m') at the bottom, decoded optimally by FHT;
**(3) Aggregation** — re-estimate the LLR of every position as the
average over all directions of `±LLR(x ⊕ z)`, sign taken from the
decoded projection's verdict. Iterations stop early once hard
decisions stabilize (`iters = 0` uses the paper's default `⌈m/2⌉`).
Complexity O(iters · n² log n) per recursion level. Message
convention: `:monomial`.

### Belief propagation (sum-product) decoder
**File:** [`src/decoders/bp.jl`](src/decoders/bp.jl) — `BPDecoder`

Sum-product decoding on the Tanner graph whose check nodes are all
rows of the dual code's generator matrix — for RM(r, m) the dual is
RM(m-r-1, m), giving a redundant (overcomplete) parity-check set,
which helps BP on these dense codes. Flooding schedule with the exact
tanh check-node rule, early exit once hard decisions satisfy every
check. BP is known to be a weak decoder for Reed-Muller codes (their
Tanner graphs have many short cycles); it is included as a comparison
baseline and as a constituent for ensemble schemes such as
`AutomorphismEnsembleDecoder`, which recovers the multiple-bases BP
idea. Message convention: `:monomial`.

### Global-list-with-permutations (GLP) decoder
**File:** [`src/decoders/glp.jl`](src/decoders/glp.jl) — `GLPDecoder`, `glp_permutations`
**References:** K. Shabunov, ecclab (`dtrm_glp` codec),
https://github.com/kshabunov/ecclab; the technique appears in Dumer &
Shabunov, "Soft-decision decoding of Reed-Muller codes: recursive
lists", 2006, §VI.

Recursive list decoding over an ensemble of code-preserving coordinate
permutations sharing **one global list**. Where
`AutomorphismEnsembleDecoder` runs its constituent decoder
independently per permutation, GLP seeds the list of a single
`DumerShabunovDecoder`-style recursion with one path per permuted
version of the received word; paths descended from different
permutations then compete for the same `L` list slots, so the budget
flows to whichever views of the channel output look most promising.
`perms` selects the permutation ensemble: `:identity` (`P1`,
equivalent to plain Dumer-Shabunov), `:cyclic` (`Pcyclic`, the m
cyclic shifts), `:pairs` (`PmCr`, every variable pair moved behind the
others), or an explicit list of permutations. Cost is one list
decoding at list size ≥ number of permutations. Message convention:
`:plotkin`.

### Local graph-search decoder
**File:** [`src/decoders/graph_search.jl`](src/decoders/graph_search.jl) — `GraphSearchDecoder`
**Reference:** M. Kamenev, "On decoding of Reed-Muller codes using a
local graph search", IEEE Trans. Commun., 2022.

Views codewords as nodes of a graph whose edges join pairs at Hamming
distance `2^(m-r)` (the minimum distance); a full traversal maximizing
the correlation `M = Σᵢ (1-2cᵢ) yᵢ` is exactly ML decoding (the
paper's Proposition 1). The decoder starts from `DumerDecoder(leaves =
:fht)`'s output and takes up to `iters` greedy steps: each step moves
to the best adjacent, not-yet-visited codeword. Adjacent codewords are
found by a greedy depth-first search down the tree of nested affine
subspaces of EG(m, 2) whose leaves enumerate all minimum-weight
codeword supports (the paper's `NextStepGreedy`): at each tree node
the child minimizing the flip penalty is chosen, with all children of
a node scored at once via an FHT, so one step costs O(n log n). The
`l` best root children are tried per step; when none yields a new
codeword, up to `lbar` extra children are tried, at most `s` times per
decode (paper defaults: l = 8, l̄ = 8, s = 5). The paper's optional
CRC-based early termination is not implemented. Message convention:
`:plotkin`.

## Generic decoder wrappers

These compose with **any** decoder — of either message basis, via the
`inner`/`code` arguments — rather than implementing a decoding
algorithm of their own.

### Automorphism Ensemble Decoding (AED)
**File:** [`src/decoders/aed.jl`](src/decoders/aed.jl) — `AutomorphismEnsembleDecoder`
**Reference:** M. Geiselhart, A. Elkelesh, M. Ebada, S. Cammerer,
S. ten Brink, "Automorphism ensemble decoding of Reed-Muller codes",
IEEE Trans. Commun., 2021.

At construction, `size` permutations of the coordinates are drawn from
the code's automorphism group (the general affine group GA(m):
x ↦ Ax ⊕ b with A invertible), the identity always included. Decoding
runs the inner decoder on each permuted version of the channel LLRs,
maps every candidate codeword back, and keeps the one with the best
correlation to the channel output — so the ensemble is never worse
than the inner decoder alone. Cost is `size` × the constituent
decoder.

### Chase-II decoding
**File:** [`src/decoders/wrappers.jl`](src/decoders/wrappers.jl) — `ChaseDecoder`
**Reference:** D. Chase, "A class of algorithms for decoding block
codes with channel measurement information", IEEE Trans. Inf. Theory,
1972.

The `t` least reliable positions are identified and all 2^t sign-flip
test patterns of the channel LLRs are decoded by the inner decoder;
the candidate codeword with the best correlation to the unmodified
channel output wins. Turns a hard-decision decoder (e.g. `ReedDecoder`)
into a genuinely soft-input one, at cost `2^t` × inner.

### Generalized Minimum Distance (GMD) decoding
**File:** [`src/decoders/wrappers.jl`](src/decoders/wrappers.jl) — `GMDDecoder`
**Reference:** G. D. Forney, "Generalized minimum distance decoding",
IEEE Trans. Inf. Theory, 1966.

The `2j` least reliable positions are erased (LLR set to 0) for
`j = 0, 1, …, ⌊(d-1)/2⌋`, each trial decoded by the inner decoder, and
the candidate with the best correlation to the channel output wins.
Includes the no-erasure trial, so it is never worse than the inner
decoder alone.

### Brute-force maximum-likelihood decoder
**File:** [`src/decoders/wrappers.jl`](src/decoders/wrappers.jl) — `MLDecoder`

Exhaustive-correlation ML decoding over all 2^k codewords — the exact
reference against which every other decoder in this package is
verified in the test suite. Only practical for small codes; the
constructor refuses k > 24. O(2^k · n) per decoded word.

## Channels and simulation

Not decoding algorithms, but part of the shared comparison harness:

- **`BSC(p)`** and **`BIAWGN(sigma)`** — [`src/channels.jl`](src/channels.jl):
  binary symmetric and binary-input AWGN channel models; `transmit`
  turns a codeword into the LLR vector a decoder expects.
- **`simulate(enc, dec, code, channel; trials)`** — [`src/simulate.jl`](src/simulate.jl):
  Monte-Carlo harness reporting bit and block error rates, used to
  compare any two encoder/decoder pipelines under identical
  conditions.

## Candidates not yet implemented

See the README's "Adding a new algorithm" section for how to
contribute one of these, or another algorithm entirely:

- Kabatiansky-Tavernier (and Fourquet-Tavernier) deterministic list
  decoding of RM(2, m) — deferred because a faithful implementation
  needs branch-and-bound details from the original paper that weren't
  available to verify against.
- Sparse/multi-decoder RPA and the hardware-oriented IPA variant.
