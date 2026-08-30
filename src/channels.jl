# Channel models.

"""
    transmit(rng::AbstractRNG, ch, codeword::BitVector) -> Vector{Float64}

Pass `codeword` through channel `ch` (a [`BSC`](@ref) or
[`BIAWGN`](@ref)) and return the resulting log-likelihood ratios
(positive = bit 0 more likely), ready for any [`AbstractDecoder`](@ref).
"""
function transmit end

"""
    BSC(p)

Binary symmetric channel with crossover probability `p`.
"""
struct BSC
    p::Float64
    function BSC(p::Real)
        0 <= p < 0.5 || throw(ArgumentError("BSC requires 0 <= p < 0.5, got $p"))
        new(p)
    end
end

function transmit(rng::AbstractRNG, ch::BSC, c::BitVector)
    scale = ch.p == 0 ? 1.0 : log((1 - ch.p) / ch.p)
    [(c[i] ⊻ (rand(rng) < ch.p)) ? -scale : scale for i in eachindex(c)]
end

"""
    BIAWGN(sigma)

Binary-input AWGN channel: BPSK modulation x = 1 - 2c, received
y = x + sigma·N(0,1), LLR = 2y/sigma². For code rate R,
Eb/N0 [dB] = 10·log10(1 / (2·R·sigma²)).
"""
struct BIAWGN
    sigma::Float64
    function BIAWGN(sigma::Real)
        sigma > 0 || throw(ArgumentError("BIAWGN requires sigma > 0, got $sigma"))
        new(sigma)
    end
end

"""
    BIAWGN_from_ebn0(ebn0_db, code) -> BIAWGN

Construct the channel from Eb/N0 in dB for the given code's rate.
"""
BIAWGN_from_ebn0(ebn0_db::Real, code::RMCode) =
    BIAWGN(sqrt(1 / (2 * rate(code) * 10.0^(ebn0_db / 10))))

transmit(rng::AbstractRNG, ch::BIAWGN, c::BitVector) =
    [2 * ((c[i] ? -1.0 : 1.0) + ch.sigma * randn(rng)) / ch.sigma^2 for i in eachindex(c)]
