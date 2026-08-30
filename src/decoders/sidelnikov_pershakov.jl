"""
    SidelnikovPershakovDecoder(; voting = :weighted)

Sidel'nikov-Pershakov derivative decoding for second-order codes
RM(2, m) (Sidel'nikov & Pershakov, "Decoding of Reed-Muller codes
with a large number of errors", Probl. Inf. Transm., 1992).

A codeword is a quadratic Boolean function f(x) = xᵀQx ⊕ ⟨l, x⟩ ⊕ a₀,
so its derivative in any direction b, D_b f(x) = f(x) ⊕ f(x ⊕ b) =
⟨Bb, x⟩ ⊕ f(b) ⊕ f(0) with B = Q + Qᵀ, is affine. The decoder:

 1. For every nonzero direction b, forms derivative LLRs by min-sum
    combining and ML-decodes them in RM(1, m) with the fast
    Walsh-Hadamard transform, giving an estimate w_b ≈ Bb with a
    reliability (the FHT peak magnitude).
 2. Recovers each entry B_ij by a reliability-weighted majority vote
    over all pairs (b, b ⊕ e_i), whose estimates XOR to B e_i.
 3. Peels the recovered quadratic part off the channel LLRs and
    FHT-decodes the residual RM(1, m) for the affine part.

`voting` selects the vote combination in step 2: `:weighted` scales
each vote by the reliability of its less reliable member, `:majority`
counts every vote equally — the simplified plain-majority voting of
Sakkour's variant of the algorithm (Sakkour, "Decoding of second
order Reed-Muller codes with a large number of errors", ITW 2005).

Corrects any pattern of up to ⌊(d-1)/2⌋ = 2^(m-3) - 1 errors and, on
random errors, most patterns well beyond half the minimum distance.
Complexity O(n² log n). Message convention: `:monomial`.
"""
struct SidelnikovPershakovDecoder <: AbstractDecoder
    voting::Symbol
    function SidelnikovPershakovDecoder(; voting::Symbol = :weighted)
        voting in (:weighted, :majority) ||
            throw(ArgumentError("voting must be :weighted or :majority, got $voting"))
        new(voting)
    end
end

basis(::SidelnikovPershakovDecoder) = :monomial

function decode(dec::SidelnikovPershakovDecoder, code::RMCode, llr::AbstractVector{<:Real})
    code.r == 2 || throw(ArgumentError("SidelnikovPershakovDecoder only handles RM(2, m), got $code"))
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    m = code.m
    L = Vector{Float64}(llr)

    # Step 1: decode the derivative in every nonzero direction.
    linmask = zeros(Int, n)      # linmask[b+1]: estimate of Bb as a bit mask
    peak = zeros(Float64, n)     # reliability of that estimate
    deriv = Vector{Float64}(undef, n)
    for b in 1:(n - 1)
        for x in 0:(n - 1)
            deriv[x + 1] = _combine(:minsum, L[x + 1], L[(x ⊻ b) + 1])
        end
        w = copy(deriv)
        _fwht!(w)
        best = argmax(abs.(w))
        linmask[b + 1] = best - 1
        peak[b + 1] = abs(w[best])
    end
    peak[1] = Inf                # D_0 f = 0 is known exactly

    # Step 2: weighted majority vote for the entries of B. Each pair
    # (b, b ⊕ e_i) votes on column B e_i with the weight of its less
    # reliable member; votes[i, j] and votes[j, i] both score B_ij.
    votes = zeros(m, m)
    for i in 1:m
        ei = 1 << (i - 1)
        for b in 0:(n - 1)
            b & ei == 0 || continue
            est = linmask[b + 1] ⊻ linmask[(b ⊻ ei) + 1]
            w = dec.voting === :weighted ? min(peak[b + 1], peak[(b ⊻ ei) + 1]) : 1.0
            for j in 1:m
                j == i && continue
                votes[i, j] += (est >> (j - 1)) & 1 == 1 ? w : -w
            end
        end
    end

    # Step 3: peel the quadratic part off and decode the affine rest.
    quadmask = falses(n)
    B = falses(m, m)
    for i in 1:m, j in (i + 1):m
        if votes[i, j] + votes[j, i] > 0
            B[i, j] = true
            quadmask .⊻= monomial_row([i, j], m)
        end
    end
    residual = [quadmask[x + 1] ? -L[x + 1] : L[x + 1] for x in 0:(n - 1)]
    _fwht!(residual)
    best = argmax(abs.(residual))
    a = best - 1

    # Assemble the monomial-basis message: constant, linear, pairs.
    msg = falses(dimension(code))
    msg[1] = residual[best] < 0
    for i in 1:m
        msg[1 + i] = (a >> (i - 1)) & 1 == 1
    end
    p = 1 + m
    for i in 1:m, j in (i + 1):m
        p += 1
        msg[p] = B[i, j]
    end
    msg
end
