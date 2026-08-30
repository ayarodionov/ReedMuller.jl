"""
    GLPDecoder(code, L; perms = :pairs, combine = :minsum, leaves = :fht)

Global-list-with-permutations decoding — recursive list decoding over
an ensemble of code-preserving coordinate permutations sharing one
global list, after the `dtrm_glp` codec of Shabunov's ecclab
(github.com/kshabunov/ecclab; the technique appears in Dumer &
Shabunov, "Soft-decision decoding of Reed-Muller codes: recursive
lists", 2006, §VI).

Where [`AutomorphismEnsembleDecoder`](@ref) runs its constituent
decoder independently per permutation, GLP seeds the list of a single
[`DumerShabunovDecoder`](@ref)-style recursion with one path per
permuted version of the received word; from then on paths descended
from different permutations compete for the same `L` list slots, so
the budget flows to whichever views of the channel output look most
promising.

`perms` selects the ensemble of permutations of the m variables
(each preserves RM(r, m)); the ecclab parameter sets are provided as
symbols, or pass an explicit `Vector` of m-element permutation
vectors:

  * `:identity` — no permutation (`P1`); equivalent to plain
    Dumer-Shabunov decoding.
  * `:cyclic`   — the m cyclic shifts of the variables (`Pcyclic`).
  * `:pairs`    — for every variable pair (i, j), the permutation
    moving i, j behind the others (`PmCr`, C(m,2) permutations).

`combine` and `leaves` are as in [`DumerShabunovDecoder`](@ref).
Cost is one list decoding at list size ≥ number of permutations.
Message convention: `:plotkin`.
"""
struct GLPDecoder <: AbstractDecoder
    inner::DumerShabunovDecoder
    pos::Vector{Vector{Int}}    # position permutations, applied to LLRs
    ipos::Vector{Vector{Int}}   # their inverses
end

function GLPDecoder(code::RMCode, L::Integer; perms = :pairs,
                    combine::Symbol = :minsum, leaves::Symbol = :fht)
    vperms = perms isa Symbol ? glp_permutations(code.m, perms) :
             [collect(Int, p) for p in perms]
    isempty(vperms) && throw(ArgumentError("permutation set is empty"))
    for p in vperms
        sort(p) == collect(1:code.m) ||
            throw(ArgumentError("$p is not a permutation of 1:$(code.m)"))
    end
    L >= length(vperms) ||
        throw(ArgumentError("list size $L is smaller than the ensemble ($(length(vperms)) permutations)"))
    pos = [_position_perm(p, code.m) for p in vperms]
    GLPDecoder(DumerShabunovDecoder(L; combine, leaves), pos, [invperm(σ) for σ in pos])
end

basis(::GLPDecoder) = :plotkin

"""
    glp_permutations(m, kind) -> Vector{Vector{Int}}

The ecclab permutation ensembles for [`GLPDecoder`](@ref): `:identity`,
`:cyclic` (m shifts) or `:pairs` (each of the C(m,2) variable pairs
moved behind the others, identity included as the last pair).
"""
function glp_permutations(m::Integer, kind::Symbol)
    if kind === :identity
        return [collect(1:m)]
    elseif kind === :cyclic
        return [circshift(collect(1:m), s) for s in 0:(m - 1)]
    elseif kind === :pairs
        return [vcat([k for k in 1:m if k != i && k != j], [i, j])
                for i in 1:m for j in (i + 1):m]
    else
        throw(ArgumentError("unknown permutation set $kind (expected :identity, :cyclic or :pairs)"))
    end
end

# The position permutation induced by a variable permutation p:
# bit i of the source index is read from bit p[i].
function _position_perm(p::Vector{Int}, m::Int)
    n = 1 << m
    σ = Vector{Int}(undef, n)
    for z in 0:(n - 1)
        y = 0
        for i in 1:m
            (z >> (p[i] - 1)) & 1 == 1 && (y |= 1 << (i - 1))
        end
        σ[z + 1] = y + 1
    end
    σ
end

function decode(dec::GLPDecoder, code::RMCode, llr::AbstractVector{<:Real})
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    L = Vector{Float64}(llr)
    paths = [_LSPath(0.0, L[σ], t) for (t, σ) in enumerate(dec.pos)]
    results = _ls_node(paths, code.r, code.m, dec.inner)
    best = results[argmin([p.metric for p in results])]
    # Map the winning codeword back through its permutation and
    # recover the message from the unpermuted (noiseless) codeword.
    cw = best.cw[dec.ipos[best.tag]]
    decode(DumerDecoder(), code, hard_llr(cw))
end
