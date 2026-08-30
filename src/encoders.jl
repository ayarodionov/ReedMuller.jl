"""
    MatrixEncoder(code; basis=:monomial)

Reference encoder: multiplies the message by a precomputed generator
matrix over F₂. Works with either message basis.
"""
struct MatrixEncoder <: AbstractEncoder
    G::BitMatrix
    basis::Symbol
end

MatrixEncoder(code::RMCode; basis::Symbol = :monomial) =
    MatrixEncoder(generator_matrix(code; basis), basis)

basis(enc::MatrixEncoder) = enc.basis

function encode(enc::MatrixEncoder, code::RMCode, message::AbstractVector)
    k, n = size(enc.G)
    length(message) == k || throw(DimensionMismatch("expected $k message bits, got $(length(message))"))
    c = falses(n)
    for j in 1:k
        if Bool(message[j])
            c .⊻= @view enc.G[j, :]
        end
    end
    c
end

"""
    PlotkinEncoder()

Fast recursive (u, u+v) encoder: O(n log n) XORs, no stored matrix.
Message convention is `:plotkin` (u-part bits first, then v-part),
identical to `MatrixEncoder(code; basis=:plotkin)`.
"""
struct PlotkinEncoder <: AbstractEncoder end

basis(::PlotkinEncoder) = :plotkin

function encode(::PlotkinEncoder, code::RMCode, message::AbstractVector)
    k = dimension(code)
    length(message) == k || throw(DimensionMismatch("expected $k message bits, got $(length(message))"))
    _plotkin_encode(BitVector(Bool.(message)), code.r, code.m)
end

function _plotkin_encode(msg::BitVector, r::Int, m::Int)
    m == 0 && return copy(msg)
    r == 0 && return msg[1] ? trues(1 << m) : falses(1 << m)
    ru = min(r, m - 1)
    ku = dimension(ru, m - 1)
    u = _plotkin_encode(msg[1:ku], ru, m - 1)
    v = _plotkin_encode(msg[(ku + 1):end], r - 1, m - 1)
    vcat(u, u .⊻ v)
end
