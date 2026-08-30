"""
    ReedDecoder()

Reed's classic majority-logic decoder (Reed, 1954). Hard-decision:
LLRs are thresholded to bits, then monomial coefficients are recovered
degree by degree, from degree r down to 0. The coefficient of each
degree-d monomial x_S is the majority over 2^(m-d) parity checks, one
per assignment of the variables outside S; after a degree is finished
its contribution is subtracted from the received word.

Guaranteed to correct up to ⌊(2^(m-r) - 1) / 2⌋ errors.
Message convention: `:monomial`.
"""
struct ReedDecoder <: AbstractDecoder end

basis(::ReedDecoder) = :monomial

function decode(::ReedDecoder, code::RMCode, llr::AbstractVector{<:Real})
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    r, m = code.r, code.m
    y = BitVector(llr .< 0)          # residual received word
    mons = monomials(r, m)
    msg = falses(length(mons))

    for d in r:-1:0
        idxs = [j for j in eachindex(mons) if length(mons[j]) == d]
        # All degree-d coefficients are voted on the same residual...
        for j in idxs
            msg[j] = _majority_vote(y, mons[j], m)
        end
        # ...then their contribution is peeled off together.
        for j in idxs
            msg[j] && (y .⊻= monomial_row(mons[j], m))
        end
    end
    msg
end

# Majority vote for the coefficient of ∏_{i∈S} x_i: one check per
# assignment b of the m-d variables outside S; each check XORs the
# residual over the 2^d points of the subcube where the S-variables
# range freely. Ties (only possible with an even number of checks)
# resolve to 0.
function _majority_vote(y::BitVector, S::Vector{Int}, m::Int)
    d = length(S)
    T = setdiff(1:m, S)
    nchecks = 1 << (m - d)
    ones_count = 0
    for b in 0:(nchecks - 1)
        base = 0
        for (t, i) in enumerate(T)
            if (b >> (t - 1)) & 1 == 1
                base |= 1 << (i - 1)
            end
        end
        check = false
        for a in 0:((1 << d) - 1)
            z = base
            for (s, i) in enumerate(S)
                if (a >> (s - 1)) & 1 == 1
                    z |= 1 << (i - 1)
                end
            end
            check ⊻= y[z + 1]
        end
        ones_count += check
    end
    2 * ones_count > nchecks
end
