"""
    ReedMuller

A common framework for implementing and comparing Reed-Muller
encoding and decoding algorithms in Julia.

Every algorithm plugs into the same interface:

  * `RMCode(r, m)` describes the code RM(r, m).
  * Encoders subtype `AbstractEncoder` and implement
    `encode(enc, code, message::BitVector)::BitVector`.
  * Decoders subtype `AbstractDecoder` and implement
    `decode(dec, code, llr::Vector{Float64})::BitVector`
    (log-likelihood ratios, positive = bit 0 more likely) and
    `basis(dec)::Symbol` naming the message-coordinate convention
    they use (`:monomial` or `:plotkin`).
  * `simulate(...)` runs encoder/channel/decoder pipelines and
    reports bit and word error rates, so different algorithms can
    be compared under identical conditions.
"""
module ReedMuller

using Random
using Printf

export RMCode, blocklength, dimension, minimum_distance, rate
export AbstractEncoder, AbstractDecoder, encode, decode, basis
export generator_matrix, monomials
export MatrixEncoder, PlotkinEncoder
export ReedDecoder, FHTDecoder, DumerDecoder, DumerShabunovDecoder
export SidelnikovPershakovDecoder, RPADecoder, BPDecoder
export AutomorphismEnsembleDecoder, ChaseDecoder, GMDDecoder, MLDecoder
export BSC, BIAWGN, BIAWGN_from_ebn0, transmit, hard_llr
export simulate, SimResult

include("code.jl")
include("matrices.jl")
include("encoders.jl")
include("decoders/reed.jl")
include("decoders/fht.jl")
include("decoders/dumer.jl")
include("decoders/dumer_shabunov.jl")
include("decoders/sidelnikov_pershakov.jl")
include("decoders/rpa.jl")
include("decoders/bp.jl")
include("decoders/aed.jl")
include("decoders/wrappers.jl")
include("channels.jl")
include("simulate.jl")

end # module
