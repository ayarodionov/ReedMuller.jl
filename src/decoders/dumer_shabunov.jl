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

With `leaves = :fht` the recursion instead terminates at first-order
nodes, where every path is expanded into the full list of affine
codewords scored by their exact correlation with the node LLRs (the
list-FHT termination of the original paper). `leaves = :bits`
(default) recurses down to repetition/single-bit leaves.

`L = 1` with `:bits` leaves reduces to plain recursive (Dumer)
decoding; growing `L` approaches maximum-likelihood performance at
cost O(L·n·log n). Message convention: `:plotkin`.
"""
struct DumerShabunovDecoder <: AbstractDecoder
    L::Int
    combine::Symbol
    leaves::Symbol
    function DumerShabunovDecoder(L::Integer = 8; combine::Symbol = :minsum,
                                  leaves::Symbol = :bits)
        L >= 1 || throw(ArgumentError("list size must be >= 1, got $L"))
        combine in (:minsum, :exact) ||
            throw(ArgumentError("combine must be :minsum or :exact, got $combine"))
        leaves in (:bits, :fht) ||
            throw(ArgumentError("leaves must be :bits or :fht, got $leaves"))
        new(L, combine, leaves)
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
    results = _ls_node([start], code.r, code.m, dec)
    results[argmin([p.metric for p in results])].msg
end

function _ls_node(paths::Vector{_LSPath}, r::Int, m::Int, dec::DumerShabunovDecoder)
    L = dec.L
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

    # List-FHT leaf: expand each path into every affine codeword of
    # RM(1, m), with the exact penalty (Σ|LLR| - correlation)/2 read
    # off the path's Hadamard spectrum; keep the L best overall.
    # (The r = 1 :plotkin message equals the monomial coefficients.)
    if r == 1 && dec.leaves === :fht
        cands = Tuple{Float64, Int, Bool, Int}[]     # metric, j, a0, tag
        for p in paths
            s1 = sum(abs, p.llr)
            w = copy(p.llr)
            _fwht!(w)
            for j in 0:(length(w) - 1)
                push!(cands, (p.metric + (s1 - w[j + 1]) / 2, j, false, p.tag))
                push!(cands, (p.metric + (s1 + w[j + 1]) / 2, j, true, p.tag))
            end
        end
        keep = length(cands) <= L ? cands : partialsort(cands, 1:L; by = first)
        return [_LSResult(met, _affine_msg(j, a0, m), _affine_cw(j, a0, m), tag)
                for (met, j, a0, tag) in keep]
    end

    half = 1 << (m - 1)

    # v phase: every path descends into RM(r-1, m-1) with combined LLRs.
    vpaths = Vector{_LSPath}(undef, length(paths))
    for (i, p) in enumerate(paths)
        Lv = [_combine(dec.combine, p.llr[j], p.llr[half + j]) for j in 1:half]
        vpaths[i] = _LSPath(p.metric, Lv, i)
    end
    vres = _ls_node(vpaths, r - 1, m - 1, dec)

    # u phase: each surviving v-hypothesis descends into RM(min(r, m-1), m-1)
    # with its v-corrected LLRs.
    upaths = Vector{_LSPath}(undef, length(vres))
    for (j, v) in enumerate(vres)
        pl = paths[v.tag].llr
        Lu = [pl[i] + (v.cw[i] ? -pl[half + i] : pl[half + i]) for i in 1:half]
        upaths[j] = _LSPath(v.metric, Lu, j)
    end
    ures = _ls_node(upaths, min(r, m - 1), m - 1, dec)

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
