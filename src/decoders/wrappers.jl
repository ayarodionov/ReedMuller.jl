# Generic algorithm-agnostic wrappers: Chase-II, GMD and brute-force
# ML. Like AutomorphismEnsembleDecoder they compose with any decoder
# through the common interface.

"""
    ChaseDecoder(code, inner; t = 4)

Chase-II decoding (Chase, "A class of algorithms for decoding block
codes with channel measurement information", 1972) around any inner
decoder. The `t` least reliable positions are identified and all 2^t
sign-flip test patterns of the channel LLRs are decoded by `inner`;
the candidate codeword with the best correlation to the unmodified
channel output wins. Turns a hard-decision decoder (e.g.
[`ReedDecoder`](@ref)) into a genuinely soft-input one at cost
2^t × inner. Message convention: that of `inner`.
"""
struct ChaseDecoder <: AbstractDecoder
    inner::AbstractDecoder
    enc::MatrixEncoder
    t::Int
end

function ChaseDecoder(code::RMCode, inner::AbstractDecoder; t::Integer = 4)
    0 <= t <= 20 || throw(ArgumentError("t must be in 0:20, got $t"))
    ChaseDecoder(inner, MatrixEncoder(code; basis = basis(inner)), t)
end

basis(dec::ChaseDecoder) = basis(dec.inner)

function decode(dec::ChaseDecoder, code::RMCode, llr::AbstractVector{<:Real})
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    L = Vector{Float64}(llr)
    weak = partialsortperm(abs.(L), 1:min(dec.t, n))
    _best_of_trials(dec.inner, dec.enc, code, L, 1 << length(weak)) do trial, L2
        for (b, i) in enumerate(weak)
            (trial - 1) >> (b - 1) & 1 == 1 && (L2[i] = -L2[i])
        end
    end
end

"""
    GMDDecoder(code, inner)

Generalized Minimum Distance decoding (Forney, 1966) around any
soft-input inner decoder: the 2j least reliable positions are erased
(LLR set to 0) for j = 0, 1, …, ⌊(d-1)/2⌋, each trial is decoded by
`inner`, and the candidate with the best correlation to the channel
output wins. Includes the no-erasure trial, so it is never worse than
`inner` alone. Message convention: that of `inner`.
"""
struct GMDDecoder <: AbstractDecoder
    inner::AbstractDecoder
    enc::MatrixEncoder
end

GMDDecoder(code::RMCode, inner::AbstractDecoder) =
    GMDDecoder(inner, MatrixEncoder(code; basis = basis(inner)))

basis(dec::GMDDecoder) = basis(dec.inner)

function decode(dec::GMDDecoder, code::RMCode, llr::AbstractVector{<:Real})
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    L = Vector{Float64}(llr)
    order = sortperm(abs.(L))
    trials = (minimum_distance(code) - 1) ÷ 2 + 1
    _best_of_trials(dec.inner, dec.enc, code, L, trials) do trial, L2
        L2[order[1:2 * (trial - 1)]] .= 0.0
    end
end

# Run `inner` on `trials` modified copies of the channel LLRs (the
# closure prepares copy `L2` for the given 1-based trial index) and
# return the message whose codeword correlates best with the channel.
function _best_of_trials(prepare!, inner::AbstractDecoder, enc::MatrixEncoder,
                         code::RMCode, L::Vector{Float64}, trials::Int)
    n = length(L)
    best_corr = -Inf
    best_msg = falses(dimension(code))
    for trial in 1:trials
        L2 = copy(L)
        prepare!(trial, L2)
        msg = decode(inner, code, L2)
        cw = encode(enc, code, msg)
        corr = sum(L[i] * (cw[i] ? -1.0 : 1.0) for i in 1:n)
        if corr > best_corr
            best_corr, best_msg = corr, msg
        end
    end
    best_msg
end

"""
    MLDecoder(code; basis = :monomial)

Brute-force maximum-likelihood decoding by exhaustive correlation
over all 2^k codewords — the exact reference against which every
other decoder can be judged. Only practical for small codes; the
constructor refuses k > 24. O(2^k · n) per decoded word.
"""
struct MLDecoder <: AbstractDecoder
    enc::MatrixEncoder
    k::Int
end

function MLDecoder(code::RMCode; basis::Symbol = :monomial)
    k = dimension(code)
    k <= 24 || throw(ArgumentError("k = $k is too large for brute-force ML (max 24)"))
    MLDecoder(MatrixEncoder(code; basis), k)
end

basis(dec::MLDecoder) = dec.enc.basis

function decode(dec::MLDecoder, code::RMCode, llr::AbstractVector{<:Real})
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    L = Vector{Float64}(llr)
    best_corr = -Inf
    best = 0
    msg = falses(dec.k)
    for w in 0:(1 << dec.k) - 1
        for i in 1:dec.k
            msg[i] = (w >> (i - 1)) & 1 == 1
        end
        cw = encode(dec.enc, code, msg)
        corr = sum(L[i] * (cw[i] ? -1.0 : 1.0) for i in 1:n)
        if corr > best_corr
            best_corr, best = corr, w
        end
    end
    for i in 1:dec.k
        msg[i] = (best >> (i - 1)) & 1 == 1
    end
    msg
end
