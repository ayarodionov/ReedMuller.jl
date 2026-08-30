"""
    RMCode(r, m)

The Reed-Muller code RM(r, m): length `2^m`, dimension
`sum(binomial(m, i) for i in 0:r)`, minimum distance `2^(m-r)`.
"""
struct RMCode
    r::Int
    m::Int
    function RMCode(r::Integer, m::Integer)
        0 <= r <= m || throw(ArgumentError("RM(r, m) requires 0 <= r <= m, got r=$r, m=$m"))
        new(r, m)
    end
end

blocklength(c::RMCode) = 1 << c.m
dimension(r::Integer, m::Integer) = sum(binomial(m, i) for i in 0:r)
dimension(c::RMCode) = dimension(c.r, c.m)
minimum_distance(c::RMCode) = 1 << (c.m - c.r)
rate(c::RMCode) = dimension(c) / blocklength(c)

Base.show(io::IO, c::RMCode) =
    print(io, "RM($(c.r), $(c.m)) [n=$(blocklength(c)), k=$(dimension(c)), d=$(minimum_distance(c))]")

"""
Abstract supertype for encoders. Implement
`encode(enc, code, message::BitVector)::BitVector` and `basis(enc)`.
"""
abstract type AbstractEncoder end

"""
Abstract supertype for decoders. Implement
`decode(dec, code, llr::Vector{Float64})::BitVector` and `basis(dec)`.

The LLR convention is `llr[i] = log(P(cᵢ=0 | yᵢ) / P(cᵢ=1 | yᵢ))`,
i.e. positive values favour bit 0. Hard-decision decoders threshold
internally (`llr .< 0`). Use [`hard_llr`](@ref) to feed a plain bit
vector to any decoder.
"""
abstract type AbstractDecoder end

function encode end
function decode end

"""
    basis(x) -> Symbol

Message-coordinate convention of an encoder or decoder: `:monomial`
(coefficients of Boolean monomials, degree-ascending) or `:plotkin`
(recursive (u, u+v) layout, u-part first). An encoder and decoder can
only be compared in a pipeline when their bases match.
"""
function basis end

"""
    hard_llr(bits) -> Vector{Float64}

Turn a hard-decision bit vector into pseudo-LLRs (+1.0 for bit 0,
-1.0 for bit 1) so hard inputs can be fed to the common
`decode(dec, code, llr)` interface.
"""
hard_llr(bits::AbstractVector) = [b ? -1.0 : 1.0 for b in Bool.(bits)]
