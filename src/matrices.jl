# Generator matrices in the two message-coordinate conventions.
#
# Coordinate convention for evaluation points: codeword position z
# (0-based, z in 0:2^m-1) corresponds to the point x ∈ F₂^m with
# x_i = bit (i-1) of z, i.e. x_1 is the least significant bit.

"""
    monomials(r, m) -> Vector{Vector{Int}}

All monomials of degree <= r in variables x_1..x_m, each given as the
sorted list of variable indices, ordered by degree ascending and
lexicographically within a degree. This fixes the `:monomial` message
convention: message bit j is the coefficient of `monomials(r, m)[j]`.
"""
function monomials(r::Integer, m::Integer)
    result = Vector{Vector{Int}}()
    for d in 0:r
        _combinations!(result, m, d)
    end
    result
end

function _combinations!(result::Vector{Vector{Int}}, m::Int, d::Int)
    cur = Int[]
    function rec(start)
        if length(cur) == d
            push!(result, copy(cur))
            return
        end
        for i in start:m
            push!(cur, i)
            rec(i + 1)
            pop!(cur)
        end
    end
    rec(1)
    return result
end

"""
    monomial_row(S, m) -> BitVector

Evaluation vector of the monomial ∏_{i∈S} x_i over all 2^m points.
"""
function monomial_row(S::Vector{Int}, m::Integer)
    n = 1 << m
    row = falses(n)
    for z in 0:(n - 1)
        row[z + 1] = all(i -> (z >> (i - 1)) & 1 == 1, S)
    end
    row
end

"""
    generator_matrix(code::RMCode; basis=:monomial) -> BitMatrix

Generator matrix of RM(r, m) in the given message convention.

  * `:monomial` — row j evaluates `monomials(r, m)[j]`.
  * `:plotkin`  — recursive `[G(r,m-1) G(r,m-1); 0 G(r-1,m-1)]`
    construction (u-part rows first, then v-part rows).
"""
function generator_matrix(code::RMCode; basis::Symbol = :monomial)
    if basis === :monomial
        mons = monomials(code.r, code.m)
        G = falses(length(mons), blocklength(code))
        for (j, S) in enumerate(mons)
            G[j, :] = monomial_row(S, code.m)
        end
        return G
    elseif basis === :plotkin
        return _plotkin_matrix(code.r, code.m)
    else
        throw(ArgumentError("unknown basis $basis (expected :monomial or :plotkin)"))
    end
end

function _plotkin_matrix(r::Int, m::Int)
    m == 0 && return trues(1, 1)
    r == 0 && return trues(1, 1 << m)
    Gu = _plotkin_matrix(min(r, m - 1), m - 1)
    Gv = _plotkin_matrix(r - 1, m - 1)
    ku, kv = size(Gu, 1), size(Gv, 1)
    half = 1 << (m - 1)
    G = falses(ku + kv, 1 << m)
    G[1:ku, 1:half] = Gu
    G[1:ku, (half + 1):end] = Gu
    G[(ku + 1):end, (half + 1):end] = Gv
    G
end
