"""
    AutomorphismEnsembleDecoder(code, inner; size = 8, rng = Random.default_rng())

Automorphism Ensemble Decoding (Geiselhart, Elkelesh, Ebada, Cammerer
& ten Brink, "Automorphism ensemble decoding of Reed-Muller codes",
2021) — a generic wrapper around any constituent decoder.

At construction, `size` permutations of the coordinates are drawn
from the automorphism group of the code (the general affine group
GA(m): x ↦ Ax ⊕ b with A invertible), the identity always included.
Decoding runs `inner` on each permuted version of the channel LLRs,
maps every candidate codeword back, and keeps the one with the best
correlation to the channel output — so the ensemble is never worse
than `inner` alone. Cost is `size` × the constituent decoder.

Message convention: that of `inner`.
"""
struct AutomorphismEnsembleDecoder <: AbstractDecoder
    inner::AbstractDecoder
    enc::MatrixEncoder
    perms::Vector{Vector{Int}}
    iperms::Vector{Vector{Int}}
end

function AutomorphismEnsembleDecoder(code::RMCode, inner::AbstractDecoder;
                                     size::Integer = 8,
                                     rng::AbstractRNG = Random.default_rng())
    size >= 1 || throw(ArgumentError("ensemble size must be >= 1, got $size"))
    n = blocklength(code)
    perms = [collect(1:n)]
    while length(perms) < size
        cols, b = _rand_affine(rng, code.m)
        σ = Vector{Int}(undef, n)
        for x in 0:(n - 1)
            y = b
            for i in 1:code.m
                (x >> (i - 1)) & 1 == 1 && (y ⊻= cols[i])
            end
            σ[x + 1] = y + 1
        end
        push!(perms, σ)
    end
    AutomorphismEnsembleDecoder(inner, MatrixEncoder(code; basis = basis(inner)),
                                perms, [invperm(p) for p in perms])
end

basis(dec::AutomorphismEnsembleDecoder) = basis(dec.inner)

function decode(dec::AutomorphismEnsembleDecoder, code::RMCode, llr::AbstractVector{<:Real})
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    L = Vector{Float64}(llr)
    best_corr = -Inf
    best_cw = falses(n)
    for (p, ip) in zip(dec.perms, dec.iperms)
        msg = decode(dec.inner, code, L[p])
        cw = encode(dec.enc, code, msg)[ip]
        corr = sum(L[i] * (cw[i] ? -1.0 : 1.0) for i in 1:n)
        if corr > best_corr
            best_corr, best_cw = corr, cw
        end
    end
    # The winning codeword is noiseless for the inner decoder, so this
    # recovers its message exactly in the inner basis.
    decode(dec.inner, code, hard_llr(best_cw))
end

# A uniformly random affine bijection of F₂^m: column masks of an
# invertible matrix A plus an offset b.
function _rand_affine(rng::AbstractRNG, m::Int)
    while true
        cols = [rand(rng, 0:(1 << m) - 1) for _ in 1:m]
        _rank_f2(cols) == m && return cols, rand(rng, 0:(1 << m) - 1)
    end
end

function _rank_f2(cols::Vector{Int})
    pivots = Int[]
    for c in cols
        for p in pivots
            c = min(c, c ⊻ p)
        end
        c != 0 && push!(pivots, c)
    end
    length(pivots)
end
