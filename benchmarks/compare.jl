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
    ("Chase-II(16) over Reed",     MatrixEncoder(CODE), ChaseDecoder(CODE, ReedDecoder(); t = 4)),
    ("BP, redundant dual checks",  MatrixEncoder(CODE), BPDecoder(CODE)),
    ("Sidelnikov-Pershakov",       MatrixEncoder(CODE), SidelnikovPershakovDecoder()),
    ("SP, majority (Sakkour)",     MatrixEncoder(CODE), SidelnikovPershakovDecoder(voting = :majority)),
    ("RPA",                        MatrixEncoder(CODE), RPADecoder()),
    ("Dumer recursive, min-sum",   PlotkinEncoder(),    DumerDecoder()),
    ("Dumer recursive, FHT leaves", PlotkinEncoder(),   DumerDecoder(leaves = :fht)),
    ("GMD over Dumer-FHT",         PlotkinEncoder(),    GMDDecoder(CODE, DumerDecoder(leaves = :fht))),
    ("AED-8 over Dumer-FHT",       PlotkinEncoder(),
     AutomorphismEnsembleDecoder(CODE, DumerDecoder(leaves = :fht);
                                 size = 8, rng = MersenneTwister(1))),
    ("Dumer-Shabunov list, L=16",  PlotkinEncoder(),    DumerShabunovDecoder(16)),
    ("D-S list L=16, FHT leaves",  PlotkinEncoder(),    DumerShabunovDecoder(16, leaves = :fht)),
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
