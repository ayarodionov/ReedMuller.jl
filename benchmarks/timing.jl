# Per-decode wall-clock timing for every decoder in the package, at
# RM(2,8), RM(2,9), and RM(2,10) — the same codes as RESULTS.md, so
# these numbers can be read alongside that report's error rates.
#
# Single-threaded (unlike the error-rate scripts): timing wants one
# decode's true cost, not throughput under contention. Each cell
# times TRIALS random messages at a fixed 2 dB Eb/N0 (the operating
# point where RESULTS.md's decoders start to visibly separate), after
# a warmup call to exclude JIT compilation.
#
#   julia --project=. benchmarks/timing.jl

using ReedMuller
using Random
using Printf

const TRIALS = parse(Int, get(ENV, "RM_TIMING_TRIALS", "50"))
const EBN0_DB = 2.0
const MS = (8, 9, 10)

function pipelines_for(code::RMCode, rng_seed::Int)
    pe, me = PlotkinEncoder(), MatrixEncoder(code)
    m = code.m
    rng = MersenneTwister(rng_seed)
    [
        ("Reed (hard)",              me, ReedDecoder()),
        ("Chase-II(t=4)/Reed",       me, ChaseDecoder(code, ReedDecoder(); t = 4)),
        ("BP",                       me, BPDecoder(code)),
        ("Sidelnikov-Pershakov",     me, SidelnikovPershakovDecoder()),
        ("SP majority (Sakkour)",    me, SidelnikovPershakovDecoder(voting = :majority)),
        ("RPA",                      me, RPADecoder()),
        ("Dumer min-sum",            pe, DumerDecoder()),
        ("Dumer FHT leaves",         pe, DumerDecoder(leaves = :fht)),
        ("GMD/Dumer-FHT",            pe, GMDDecoder(code, DumerDecoder(leaves = :fht))),
        ("AED-8/Dumer-FHT",          pe, AutomorphismEnsembleDecoder(code, DumerDecoder(leaves = :fht); size = 8, rng)),
        ("DS16 FHT leaves",          pe, DumerShabunovDecoder(16, leaves = :fht)),
        ("GLP cyclic",               pe, GLPDecoder(code, 4 * m; perms = :cyclic)),
        ("GLP cyclic L=256",         pe, GLPDecoder(code, 256; perms = :cyclic)),
        ("GLP pairs (PmCr)",         pe, GLPDecoder(code, 2 * binomial(m, 2); perms = :pairs)),
        ("Graph search N=32",        pe, GraphSearchDecoder(iters = 32)),
    ]
end

# Mean decode time in milliseconds over TRIALS random messages at a
# fixed Eb/N0, after one untimed warmup call.
function time_decoder(enc, dec, code::RMCode, rng::AbstractRNG)
    ch = ReedMuller.BIAWGN_from_ebn0(EBN0_DB, code)
    k = dimension(code)
    warm_llr = transmit(rng, ch, encode(enc, code, bitrand(rng, k)))
    decode(dec, code, warm_llr)   # warmup: exclude JIT compilation

    total = 0.0
    for _ in 1:TRIALS
        llr = transmit(rng, ch, encode(enc, code, bitrand(rng, k)))
        total += @elapsed decode(dec, code, llr)
    end
    1000 * total / TRIALS
end

println("Single-threaded, $TRIALS decodes/cell, Eb/N0 = $EBN0_DB dB\n")
@printf("%-24s", "Decoder (ms/decode)")
foreach(m -> @printf("%14s", "RM(2,$m)"), MS)
println()
flush(stdout)

rows = Dict{String, Vector{Float64}}()
order = String[]
for m in MS
    code = RMCode(2, m)
    rng = MersenneTwister(20_000 + m)
    for (name, enc, dec) in pipelines_for(code, 30_000 + m)
        name in keys(rows) || (rows[name] = Float64[]; push!(order, name))
        push!(rows[name], time_decoder(enc, dec, code, rng))
    end
end

for name in order
    @printf("%-24s", name)
    foreach(t -> @printf("%14.3f", t), rows[name])
    println()
    flush(stdout)
end

println("\nTIMING_DONE")
flush(stdout)
