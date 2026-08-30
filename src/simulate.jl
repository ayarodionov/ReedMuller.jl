"""
    SimResult

Outcome of a `simulate` run: message bit error rate `ber`, block
(word) error rate `wer`, and the number of `trials`.
"""
struct SimResult
    ber::Float64
    wer::Float64
    trials::Int
end

Base.show(io::IO, s::SimResult) =
    @printf(io, "BER = %.3e, WER = %.3e (%d trials)", s.ber, s.wer, s.trials)

"""
    simulate(enc, dec, code, channel; trials = 10_000, rng = Random.default_rng())

Monte-Carlo comparison harness: random message → `encode` → `transmit`
→ `decode`, counting message bit and block errors. Throws if the
encoder and decoder use different message bases (their pipelines would
not be comparable).
"""
function simulate(enc::AbstractEncoder, dec::AbstractDecoder, code::RMCode, channel;
                  trials::Integer = 10_000, rng::AbstractRNG = Random.default_rng())
    basis(enc) === basis(dec) ||
        throw(ArgumentError("encoder basis $(basis(enc)) != decoder basis $(basis(dec))"))
    k = dimension(code)
    bit_errors = 0
    word_errors = 0
    for _ in 1:trials
        msg = bitrand(rng, k)
        c = encode(enc, code, msg)
        llr = transmit(rng, channel, c)
        est = decode(dec, code, llr)
        e = count(msg .⊻ est)
        bit_errors += e
        word_errors += e > 0
    end
    SimResult(bit_errors / (trials * k), word_errors / trials, trials)
end
