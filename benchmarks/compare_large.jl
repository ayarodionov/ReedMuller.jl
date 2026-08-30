# Large-scale comparison: RM(2,8), RM(2,9), RM(2,10) over BI-AWGN,
# multi-threaded (each trial runs on its own thread-local RNG).
#
# Run with multiple threads to make 100_000 trials/point tractable:
#   julia --project=. -t auto benchmarks/compare_large.jl
#
# Eb/N0 points are restricted to the transition region (0-3 dB) for
# these low-rate codes; beyond that every near-ML decoder already
# saturates the trial budget at zero errors.

using ReedMuller
using Random
using Printf
using Base.Threads

const MS = (8, 9, 10)
const TRIALS = parse(Int, get(ENV, "RM_TRIALS", "100000"))
const EBN0_DB = (0.0, 1.0, 2.0, 3.0)

# Threaded Monte-Carlo harness: statically partitions trials across
# threads, each with its own thread-local RNG (no shared mutable
# state is touched inside decode/encode, so this is safe).
function psimulate(enc, dec, code, ch; trials::Int)
    k = dimension(code)
    nt = Threads.maxthreadid()   # threadid() can exceed nthreads() (interactive/GC pools)
    bit_err = zeros(Int, nt)
    word_err = zeros(Int, nt)
    Threads.@threads :static for i in 1:trials
        tid = threadid()
        rng = Random.default_rng(tid)
        msg = bitrand(rng, k)
        c = encode(enc, code, msg)
        llr = transmit(rng, ch, c)
        est = decode(dec, code, llr)
        e = count(msg .⊻ est)
        bit_err[tid] += e
        word_err[tid] += e > 0
    end
    SimResult(sum(bit_err) / (trials * k), sum(word_err) / trials, trials)
end

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
        ("GLP pairs (PmCr)",         pe, GLPDecoder(code, 2 * binomial(m, 2); perms = :pairs)),
    ]
end

println("Threads: $(nthreads()), trials/point: $TRIALS, Eb/N0 points: $(collect(EBN0_DB))")
flush(stdout)

for m in MS
    code = RMCode(2, m)
    println("\n=== RM(2, $m): n=$(blocklength(code)), k=$(dimension(code)), d=$(minimum_distance(code)) ===")
    @printf("%-24s", "Eb/N0 [dB]")
    foreach(e -> @printf("%12.1f", e), EBN0_DB)
    println()
    flush(stdout)

    for (name, enc, dec) in pipelines_for(code, 10_000 + m)
        @printf("%-24s", name)
        flush(stdout)
        for ebn0 in EBN0_DB
            ch = ReedMuller.BIAWGN_from_ebn0(ebn0, code)
            res = psimulate(enc, dec, code, ch; trials = TRIALS)
            @printf("%12.2e", res.wer)
            flush(stdout)
        end
        println()
        flush(stdout)
    end
end

println("\nDONE")
flush(stdout)
