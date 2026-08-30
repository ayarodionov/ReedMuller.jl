"""
    DumerDecoder(; combine = :minsum)

Dumer's recursive soft-decision decoder (Dumer, "Recursive decoding
and its performance for low-rate Reed-Muller codes", 2004). Splits
RM(r, m) via the Plotkin (u, u+v) decomposition, first decoding the
v-part from combined half-LLRs, then the u-part from the v-corrected
sum. Leaf nodes: repetition codes (r = 0) decide by the LLR sum,
rate-1 codes (r = m) by symbol-wise hard decision.

`combine` chooses the check-node rule for the v-branch:
`:minsum` (sign-min approximation, the usual practical choice) or
`:exact` (2·atanh(tanh·tanh), true LLR combining).

`leaves = :fht` stops the recursion at first-order nodes and decodes
them optimally with the fast Hadamard transform instead of recursing
down to single bits (the stronger termination Dumer recommends);
`:bits` (default) recurses all the way down.

Complexity O(n log n). Message convention: `:plotkin`.
"""
struct DumerDecoder <: AbstractDecoder
    combine::Symbol
    leaves::Symbol
    function DumerDecoder(; combine::Symbol = :minsum, leaves::Symbol = :bits)
        combine in (:minsum, :exact) ||
            throw(ArgumentError("combine must be :minsum or :exact, got $combine"))
        leaves in (:bits, :fht) ||
            throw(ArgumentError("leaves must be :bits or :fht, got $leaves"))
        new(combine, leaves)
    end
end

basis(::DumerDecoder) = :plotkin

function decode(dec::DumerDecoder, code::RMCode, llr::AbstractVector{<:Real})
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    msg, _ = _dumer(dec, Vector{Float64}(llr), code.r, code.m)
    msg
end

# Returns (message bits, decoded codeword) for the RM(r, m) node.
function _dumer(dec::DumerDecoder, llr::Vector{Float64}, r::Int, m::Int)
    if r == 0
        bit = sum(llr) < 0
        return BitVector([bit]), (bit ? trues(1 << m) : falses(1 << m))
    end
    # For r = 1 the :plotkin message of RM(1, m) coincides with the
    # monomial coefficient vector [a₀, a₁, …, a_m] (the v-bit of the
    # (u, u+v) split is exactly the coefficient of the top variable),
    # so the FHT result can be returned as is.
    if r == 1 && dec.leaves === :fht
        return _fht_ml(llr, m)
    end
    if r == m
        cw = BitVector(llr .< 0)
        return _invert_full(cw, m), cw
    end
    half = 1 << (m - 1)
    L = @view llr[1:half]
    R = @view llr[(half + 1):end]

    Lv = [_combine(dec.combine, L[i], R[i]) for i in 1:half]
    mv, v = _dumer(dec, Lv, r - 1, m - 1)

    Lu = [L[i] + (v[i] ? -R[i] : R[i]) for i in 1:half]
    mu, u = _dumer(dec, Lu, r, m - 1)

    vcat(mu, mv), vcat(u, u .⊻ v)
end

# Check-node LLR combining rule, shared with DumerShabunovDecoder.
_combine(combine::Symbol, a::Float64, b::Float64) =
    combine === :minsum ? sign(a) * sign(b) * min(abs(a), abs(b)) :
                          2 * atanh(clamp(tanh(a / 2) * tanh(b / 2), -1 + 1e-15, 1 - 1e-15))

# Rate-1 node: recover the plotkin-basis message from a full-space
# codeword of RM(m, m) by unwinding u = left, v = left ⊻ right.
function _invert_full(c::BitVector, m::Int)
    m == 0 && return copy(c)
    half = 1 << (m - 1)
    u = c[1:half]
    v = u .⊻ c[(half + 1):end]
    vcat(_invert_full(u, m - 1), _invert_full(v, m - 1))
end
