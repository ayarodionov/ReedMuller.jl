"""
    DumerShabunovDecoder(L = 8; combine = :minsum)

Dumer-Shabunov recursive list decoding (Dumer & Shabunov,
"Soft-decision decoding of Reed-Muller codes: recursive lists", 2006).

Same Plotkin (u, u+v) recursion as [`DumerDecoder`](@ref), but instead
of committing to one decision per node it maintains a list of up to
`L` candidate decoding paths. Information decisions are made at the
leaves — repetition nodes (r = 0) and single-position nodes (m = 0),
to which rate-1 subcodes reduce recursively — where every surviving
path branches into both hypotheses. Paths are ranked by the standard
successive-cancellation path metric: the accumulated magnitude of
leaf LLRs that the path's decisions contradict. After each branching
the list is pruned back to the `L` best paths; the best final path is
returned.

`L = 1` reduces to plain recursive (Dumer) decoding; growing `L`
approaches maximum-likelihood performance at cost O(L·n·log n).
Message convention: `:plotkin`.
"""
struct DumerShabunovDecoder <: AbstractDecoder
    L::Int
    combine::Symbol
    function DumerShabunovDecoder(L::Integer = 8; combine::Symbol = :minsum)
        L >= 1 || throw(ArgumentError("list size must be >= 1, got $L"))
        combine in (:minsum, :exact) ||
            throw(ArgumentError("combine must be :minsum or :exact, got $combine"))
        new(L, combine)
    end
end

basis(::DumerShabunovDecoder) = :plotkin

# A live path entering a node: accumulated metric, this node's LLRs,
# and the index of the path's context in the caller's list.
struct _LSPath
    metric::Float64
    llr::Vector{Float64}
    tag::Int
end

# A path leaving a node: extended metric, the message bits and
# codeword decoded for this node, and the caller-context tag.
struct _LSResult
    metric::Float64
    msg::BitVector
    cw::BitVector
    tag::Int
end

function decode(dec::DumerShabunovDecoder, code::RMCode, llr::AbstractVector{<:Real})
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    start = _LSPath(0.0, Vector{Float64}(llr), 1)
    results = _ls_node([start], code.r, code.m, dec.L, dec.combine)
    results[argmin([p.metric for p in results])].msg
end

function _ls_node(paths::Vector{_LSPath}, r::Int, m::Int, L::Int, combine::Symbol)
    # Leaf: repetition code (covers the single-position node m = 0).
    # One information bit; every path branches into both hypotheses,
    # paying the magnitude of each disagreeing LLR.
    if r == 0 || m == 0
        n = 1 << m
        out = Vector{_LSResult}(undef, 2 * length(paths))
        for (i, p) in enumerate(paths)
            pen0 = 0.0
            pen1 = 0.0
            for l in p.llr
                l < 0 ? (pen0 -= l) : (pen1 += l)
            end
            out[2i - 1] = _LSResult(p.metric + pen0, falses(1), falses(n), p.tag)
            out[2i]     = _LSResult(p.metric + pen1, trues(1), trues(n), p.tag)
        end
        return _ls_prune(out, L)
    end

    half = 1 << (m - 1)

    # v phase: every path descends into RM(r-1, m-1) with combined LLRs.
    vpaths = Vector{_LSPath}(undef, length(paths))
    for (i, p) in enumerate(paths)
        Lv = [_combine(combine, p.llr[j], p.llr[half + j]) for j in 1:half]
        vpaths[i] = _LSPath(p.metric, Lv, i)
    end
    vres = _ls_node(vpaths, r - 1, m - 1, L, combine)

    # u phase: each surviving v-hypothesis descends into RM(min(r, m-1), m-1)
    # with its v-corrected LLRs.
    upaths = Vector{_LSPath}(undef, length(vres))
    for (j, v) in enumerate(vres)
        pl = paths[v.tag].llr
        Lu = [pl[i] + (v.cw[i] ? -pl[half + i] : pl[half + i]) for i in 1:half]
        upaths[j] = _LSPath(v.metric, Lu, j)
    end
    ures = _ls_node(upaths, min(r, m - 1), m - 1, L, combine)

    out = Vector{_LSResult}(undef, length(ures))
    for (idx, u) in enumerate(ures)
        v = vres[u.tag]
        out[idx] = _LSResult(u.metric, vcat(u.msg, v.msg),
                             vcat(u.cw, u.cw .⊻ v.cw), v.tag)
    end
    out
end

function _ls_prune(paths::Vector{_LSResult}, L::Int)
    length(paths) <= L && return paths
    partialsort(paths, 1:L; by = p -> p.metric)
end
