"""
    GraphSearchDecoder(; iters = 64, l = 8, lbar = 8, s = 5)

Local graph-search decoding (Kamenev, "On decoding of Reed-Muller
codes using a local graph search", IEEE Trans. Commun., 2022).

The codewords of RM(r, m) form a graph whose edges connect pairs at
Hamming distance `2^(m-r)` (the minimum distance); a full traversal
maximizing the correlation metric `M = Σᵢ (1-2cᵢ) yᵢ` is exactly ML
decoding (Proposition 1 of the paper). The decoder starts from the
output of Dumer's recursive decoder and takes up to `iters` greedy
steps: each step moves to an adjacent, not-yet-visited codeword —
found by flipping the support of a minimum-weight codeword — and the
best codeword seen is returned.

Adjacent codewords are found by a greedy depth-first search down the
tree of nested affine subspaces of EG(m, 2) whose leaves enumerate
all minimum-weight-codeword supports (the paper's `NextStepGreedy`):
at each tree node the child minimizing the flip penalty
`2 Σ_{i∈S̃} (1-2cᵢ) yᵢ` is chosen, with all children scored at once
by a fast Hadamard transform, so one step costs O(n log n). The `l`
best children of the root are tried per step; when none yields a new
codeword, up to `lbar` extra children are tried, at most `s` times
per decode (paper defaults: l = 8, lbar = 8, s = 5).

The paper's optional CRC-based early termination is not implemented
(this package has no CRC layer); the search always runs until
`iters` steps, no new adjacent codeword can be found, or the extra
budget is exhausted. Message convention: `:plotkin`.
"""
struct GraphSearchDecoder <: AbstractDecoder
    iters::Int
    l::Int
    lbar::Int
    s::Int
    function GraphSearchDecoder(; iters::Integer = 64, l::Integer = 8,
                                lbar::Integer = 8, s::Integer = 5)
        iters >= 1 || throw(ArgumentError("iters must be >= 1, got $iters"))
        l >= 1 || throw(ArgumentError("l must be >= 1, got $l"))
        lbar >= 0 || throw(ArgumentError("lbar must be >= 0, got $lbar"))
        s >= 0 || throw(ArgumentError("s must be >= 0, got $s"))
        new(iters, l, lbar, s)
    end
end

basis(::GraphSearchDecoder) = :plotkin

function decode(dec::GraphSearchDecoder, code::RMCode, llr::AbstractVector{<:Real})
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    y = Vector{Float64}(llr)

    init = decode(DumerDecoder(leaves = :fht), code, y)
    code.r == 0 && return init          # 2 codewords; initial decode is already ML
    c = encode(PlotkinEncoder(), code, init)

    corr(cw) = sum(cw[i] ? -y[i] : y[i] for i in 1:n)
    visited = Set{BitVector}((copy(c),))
    best_cw = copy(c)
    best_M = corr(c)
    extra_budget = dec.s
    lbar = dec.lbar

    for _ in 1:dec.iters
        cnew, used_extra = _gs_next_step(code.r, code.m, y, c, visited, dec.l, lbar)
        if used_extra
            extra_budget -= 1
            extra_budget <= 0 && (lbar = 0)
        end
        cnew === nothing && break
        c = cnew
        push!(visited, copy(c))
        M = corr(c)
        if M > best_M
            best_M = M
            best_cw = copy(c)
        end
    end
    decode(DumerDecoder(), code, hard_llr(best_cw))
end

# One greedy step: return the adjacent codeword (Hamming distance
# 2^(m-r) from c) with the smallest flip penalty among up to l (+lbar)
# greedy depth-first searches, or nothing if every attempt landed on a
# visited codeword. Root children are the supports of the 2n-2
# weight-n/2 codewords a ⊕ ⟨t, z⟩ of RM(1, m), scored via the FHT.
function _gs_next_step(r::Int, m::Int, y::Vector{Float64}, c::BitVector,
                       visited::Set{BitVector}, l::Int, lbar::Int)
    n = 1 << m
    h = [c[z + 1] ? -y[z + 1] : y[z + 1] for z in 0:(n - 1)]
    _fwht!(h)
    cands = Vector{Tuple{Float64, Int, Bool}}(undef, 2n - 2)
    for t in 1:(n - 1)
        cands[2t - 1] = (h[1] - h[t + 1], t, false)
        cands[2t]     = (h[1] + h[t + 1], t, true)
    end
    sort!(cands; by = first)

    best_flip = Inf
    best_res = nothing
    used_extra = false
    for (idx, (val, t, a)) in enumerate(cands)
        if idx > l
            best_res === nothing || break     # found within l children: stop
            (lbar == 0 || idx > l + lbar) && break
            used_extra = true
        end
        v = [z + 1 for z in 0:(n - 1) if a ⊻ isodd(count_ones(t & z))]
        res, flip = _gs_rec(r - 1, m - 1, y, c, visited, v, val)
        if res !== nothing && flip < best_flip
            best_flip = flip
            best_res = res
        end
        best_flip < 0 && break                # strict improvement found
    end
    best_res, used_extra
end

# Greedy DFS on the shortened code living on positions v (|v| = 2^m):
# descend to the child subspace with the smallest flip penalty until
# the leaf (a minimum-weight support) is reached, then flip it.
function _gs_rec(r::Int, m::Int, y::Vector{Float64}, c::BitVector,
                 visited::Set{BitVector}, v::Vector{Int}, val::Float64)
    if r == 0
        res = copy(c)
        for p in v
            res[p] = !res[p]
        end
        res in visited && return nothing, Inf
        return res, val
    end
    nn = length(v)
    h = [c[v[j]] ? -y[v[j]] : y[v[j]] for j in 1:nn]
    _fwht!(h)
    best = Inf
    bt, ba = 1, false
    for t in 1:(nn - 1)
        for a in (false, true)
            childval = a ? h[1] + h[t + 1] : h[1] - h[t + 1]
            if childval < best
                best = childval
                bt, ba = t, a
            end
        end
    end
    vhat = [v[j] for j in 1:nn if ba ⊻ isodd(count_ones(bt & (j - 1)))]
    _gs_rec(r - 1, m - 1, y, c, visited, vhat, best)
end
