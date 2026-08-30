"""
    FHTDecoder()

Maximum-likelihood soft-decision decoder for first-order codes
RM(1, m) via the fast Walsh-Hadamard transform (the "Green machine").
Complexity O(n log n). The transform of the LLR vector is maximised in
absolute value at the index whose bits are the linear coefficients;
the sign gives the constant term.

Message convention: `:monomial`, i.e. `[a₀, a₁, …, a_m]` with
codeword c(z) = a₀ ⊕ ⊕ᵢ aᵢ zᵢ.
"""
struct FHTDecoder <: AbstractDecoder end

basis(::FHTDecoder) = :monomial

function decode(::FHTDecoder, code::RMCode, llr::AbstractVector{<:Real})
    code.r == 1 || throw(ArgumentError("FHTDecoder only handles RM(1, m), got $code"))
    n = blocklength(code)
    length(llr) == n || throw(DimensionMismatch("expected $n LLRs, got $(length(llr))"))
    _fht_ml(Vector{Float64}(llr), code.m)[1]
end

# Message [a₀, a₁, …, a_m] and codeword of the affine function
# a₀ ⊕ ⟨j, x⟩, where the bits of j are the linear coefficients.
function _affine_msg(j::Int, a0::Bool, m::Int)
    msg = falses(m + 1)
    msg[1] = a0
    for i in 1:m
        msg[1 + i] = (j >> (i - 1)) & 1 == 1
    end
    msg
end

function _affine_cw(j::Int, a0::Bool, m::Int)
    n = 1 << m
    cw = falses(n)
    for x in 0:(n - 1)
        cw[x + 1] = a0 ⊻ isodd(count_ones(x & j))
    end
    cw
end

# ML decoding of RM(1, m) from LLRs: (message, codeword).
function _fht_ml(llr::Vector{Float64}, m::Int)
    w = copy(llr)
    _fwht!(w)
    best = argmax(abs.(w))
    j = best - 1
    a0 = w[best] < 0
    _affine_msg(j, a0, m), _affine_cw(j, a0, m)
end

# In-place fast Walsh-Hadamard transform:
# w[j+1] <- Σ_z w[z+1] * (-1)^⟨z, j⟩
function _fwht!(w::Vector{Float64})
    n = length(w)
    h = 1
    while h < n
        for base in 0:(2h):(n - 1)
            for i in (base + 1):(base + h)
                a, b = w[i], w[i + h]
                w[i] = a + b
                w[i + h] = a - b
            end
        end
        h *= 2
    end
    w
end
