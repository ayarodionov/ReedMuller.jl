"""
    BPDecoder(code; iters = 30)

Belief-propagation (sum-product) decoding on the Tanner graph whose
check nodes are all rows of the dual code's generator matrix — for
RM(r, m) the dual is RM(m-r-1, m), giving a redundant (overcomplete)
parity-check set, which helps BP on these dense codes. Flooding
schedule with the exact tanh check-node rule, early exit as soon as
the hard decisions satisfy every check. The final hard estimate is
mapped to message bits with a run of Reed's decoder.

BP is known to be a weak decoder for Reed-Muller codes (their Tanner
graphs have many short cycles) — it is included as a comparison
baseline and as a constituent for ensemble schemes such as
[`AutomorphismEnsembleDecoder`](@ref), which recovers the
multiple-bases BP idea. Message convention: `:monomial`.
"""
struct BPDecoder <: AbstractDecoder
    checks::Vector{Vector{Int}}   # variable indices per check node
    iters::Int
end

function BPDecoder(code::RMCode; iters::Integer = 30)
    iters >= 1 || throw(ArgumentError("iters must be >= 1, got $iters"))
    checks = Vector{Int}[]
    if code.r < code.m               # r == m has an empty dual
        H = generator_matrix(RMCode(code.m - code.r - 1, code.m))
        checks = [findall(H[i, :]) for i in axes(H, 1)]
    end
    BPDecoder(checks, iters)
end

basis(::BPDecoder) = :monomial

function decode(dec::BPDecoder, code::RMCode, llr::AbstractVector{<:Real})
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    L = Vector{Float64}(llr)
    total = copy(L)
    msg_cv = [zeros(length(c)) for c in dec.checks]

    for _ in 1:dec.iters
        hard = BitVector(total .< 0)
        all(c -> iseven(count(@view hard[c])), dec.checks) && break
        for (ci, c) in enumerate(dec.checks)
            out = msg_cv[ci]
            # Prefix/suffix products of tanh(q/2) give each variable
            # the product over all *other* variables of the check.
            th = [tanh(clamp(total[v] - out[k], -30.0, 30.0) / 2)
                  for (k, v) in enumerate(c)]
            suffix = ones(length(c) + 1)
            for k in length(c):-1:1
                suffix[k] = suffix[k + 1] * th[k]
            end
            prefix = 1.0
            for k in eachindex(c)
                rest = prefix * suffix[k + 1]
                out[k] = 2 * atanh(clamp(rest, -1 + 1e-12, 1 - 1e-12))
                prefix *= th[k]
            end
        end
        total .= L
        for (ci, c) in enumerate(dec.checks), (k, v) in enumerate(c)
            total[v] += msg_cv[ci][k]
        end
    end
    decode(ReedDecoder(), code, hard_llr(BitVector(total .< 0)))
end
