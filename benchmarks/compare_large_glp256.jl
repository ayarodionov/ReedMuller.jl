# GLPDecoder(:cyclic, L=256) row for the RM(2,8-10) comparison in
# ../RESULTS.md — ecclab's own reported configuration for these code
# lengths (Pcyclic, L=256), run separately from the main sweep
# (which uses L=4m) since it's a much larger list.
#
#   julia --project=. -t auto benchmarks/compare_large_glp256.jl

using ReedMuller
using Random
using Printf
using Base.Threads

const TRIALS = parse(Int, get(ENV, "RM_TRIALS", "100000"))
const EBN0_DB = (0.0, 1.0, 2.0, 3.0)
const MS = (8, 9, 10)
const L = 256

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
    dec = GLPDecoder(code, L; perms = :cyclic)
    println("=== RM(2, $m): n=$(blocklength(code)) ===")
    @printf("%-24s", "GLP cyclic L=256")
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

println("\nGLP256_DONE")
flush(stdout)
