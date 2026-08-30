# GraphSearchDecoder rows for the RM(2,8-10) comparison in
# ../RESULTS.md — run separately since it was added to the package
# after the main sweep completed.
#
#   julia --project=. -t auto benchmarks/compare_large_gs.jl

using ReedMuller
using Random
using Printf
using Base.Threads

const TRIALS = parse(Int, get(ENV, "RM_TRIALS", "100000"))
const EBN0_DB = (0.0, 1.0, 2.0, 3.0)
const MS = (8, 9, 10)

function psimulate(enc, dec, code, ch; trials::Int)
    k = dimension(code)
    nt = Threads.maxthreadid()
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

println("Threads: $(nthreads()), trials/point: $TRIALS, Eb/N0 points: $(collect(EBN0_DB))")
flush(stdout)

for m in MS
    code = RMCode(2, m)
    pe = PlotkinEncoder()
    dec = GraphSearchDecoder(iters = 32)
    println("=== RM(2, $m): n=$(blocklength(code)) ===")
    @printf("%-24s", "Graph search, N=32")
    flush(stdout)
    for ebn0 in EBN0_DB
        ch = ReedMuller.BIAWGN_from_ebn0(ebn0, code)
        res = psimulate(pe, dec, code, ch; trials = TRIALS)
        @printf("%12.2e", res.wer)
        flush(stdout)
    end
    println()
    flush(stdout)
end

println("\nGS_DONE")
flush(stdout)
