"""
    RPADecoder(; iters = 0)

Recursive Projection-Aggregation decoding (Ye & Abbe, "Recursive
projection-aggregation decoding of Reed-Muller codes", 2019), the
soft (LLR) variant. Performance is close to maximum likelihood for
low orders (r ≤ 3).

For RM(r, m) with r ≥ 2, each iteration:

 1. **Projection**: for every nonzero direction z the received LLRs
    are projected onto the quotient by the subspace {0, z} — the LLR
    of y(x) ⊕ y(x ⊕ z) per coset — which carries a codeword of
    RM(r-1, m-1).
 2. **Recursion**: each projection is decoded recursively; RM(1, m')
    is reached at the bottom and decoded optimally by FHT.
 3. **Aggregation**: the LLR of every position is re-estimated as the
    average over all directions of ±LLR(x ⊕ z), the sign taken from
    the decoded projection's verdict on whether y(x) = y(x ⊕ z).

Iterations stop early once the hard decisions are stable; `iters = 0`
uses the paper's default of ⌈m/2⌉. The final hard estimate is turned
back into message bits with a run of Reed's decoder. Complexity
O(iters · n² log n) per recursion level. Message convention:
`:monomial`.
"""
struct RPADecoder <: AbstractDecoder
    iters::Int
    function RPADecoder(; iters::Integer = 0)
        iters >= 0 || throw(ArgumentError("iters must be >= 0, got $iters"))
        new(iters)
    end
end

basis(::RPADecoder) = :monomial

function decode(dec::RPADecoder, code::RMCode, llr::AbstractVector{<:Real})
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    L = Vector{Float64}(llr)
    code.r == 0 && return BitVector([sum(L) < 0])
    iters = dec.iters == 0 ? cld(code.m, 2) : dec.iters
    cw = _rpa(L, code.r, code.m, iters)
    decode(ReedDecoder(), code, hard_llr(cw))
end

# Remove / re-insert bit position h (0-based) of an index: the
# quotient by {0, z} is indexed by dropping the top set bit of z from
# the coset representative.
_delbit(x::Int, h::Int) = (x & ((1 << h) - 1)) | ((x >> (h + 1)) << h)
_insbit(q::Int, h::Int) = (q & ((1 << h) - 1)) | ((q >> h) << (h + 1))

function _rpa(L::Vector{Float64}, r::Int, m::Int, iters::Int)
    r == 1 && return _fht_ml(L, m)[2]
    n = 1 << m
    half = n >> 1
    Lcur = copy(L)
    prev = BitVector(Lcur .< 0)
    proj = Vector{Float64}(undef, half)
    for _ in 1:iters
        Lnew = zeros(n)
        for z in 1:(n - 1)
            h = 8 * sizeof(Int) - leading_zeros(z) - 1   # top set bit of z
            for q in 0:(half - 1)
                x = _insbit(q, h)                        # coset representative
                proj[q + 1] = _combine(:minsum, Lcur[x + 1], Lcur[(x ⊻ z) + 1])
            end
            yz = _rpa(copy(proj), r - 1, m - 1, iters)
            for x in 0:(n - 1)
                rep = (x >> h) & 1 == 1 ? x ⊻ z : x
                flip = yz[_delbit(rep, h) + 1]
                Lnew[x + 1] += flip ? -Lcur[(x ⊻ z) + 1] : Lcur[(x ⊻ z) + 1]
            end
        end
        Lnew ./= n - 1
        Lcur = Lnew
        hard = BitVector(Lcur .< 0)
        hard == prev && break
        prev = hard
    end
    BitVector(Lcur .< 0)
end
