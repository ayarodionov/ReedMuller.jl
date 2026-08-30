# Compare all decoder pipelines on a common code and channel sweep.
#
# Run from the package directory:
#   julia --project=. benchmarks/compare.jl

using ReedMuller
using Random
using Printf

const CODE = RMCode(2, 6)          # [64, 22, 16]
const TRIALS = 5_000
const EBN0_DB = 0.0:1.0:5.0

pipelines = [
    ("Reed majority-logic (hard)", MatrixEncoder(CODE), ReedDecoder()),
    ("Sidelnikov-Pershakov",       MatrixEncoder(CODE), SidelnikovPershakovDecoder()),
    ("Dumer recursive, min-sum",   PlotkinEncoder(),    DumerDecoder()),
    ("Dumer recursive, exact",     PlotkinEncoder(),    DumerDecoder(combine = :exact)),
    ("Dumer-Shabunov list, L=4",   PlotkinEncoder(),    DumerShabunovDecoder(4)),
    ("Dumer-Shabunov list, L=16",  PlotkinEncoder(),    DumerShabunovDecoder(16)),
]

println("Code: $CODE, BI-AWGN, $TRIALS trials per point\n")
@printf("%-28s", "Eb/N0 [dB]")
foreach(e -> @printf("%12.1f", e), EBN0_DB)
println()

for (name, enc, dec) in pipelines
    rng = MersenneTwister(2026)
    @printf("%-28s", name)
    for ebn0 in EBN0_DB
        ch = ReedMuller.BIAWGN_from_ebn0(ebn0, CODE)
        res = simulate(enc, dec, CODE, ch; trials = TRIALS, rng)
        @printf("%12.2e", res.wer)
    end
    println()
end
println("\n(values are word error rates)")
