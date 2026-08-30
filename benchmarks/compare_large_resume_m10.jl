# Resume of compare_large.jl for RM(2,10) only, picking up after the
# original run died (OOM) partway through the RPA row. Reed, Chase-II,
# BP, Sidelnikov-Pershakov and SP-majority for RM(2,10) already
# completed successfully in benchmarks/results/compare_large_m8-10.log;
# this script re-runs only the remaining rows: RPA onward.
#
#   julia --project=. -t auto benchmarks/compare_large_resume_m10.jl

using ReedMuller
using Random
using Printf
using Base.Threads

const TRIALS = parse(Int, get(ENV, "RM_TRIALS", "100000"))
const EBN0_DB = (0.0, 1.0, 2.0, 3.0)

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

code = RMCode(2, 10)
m = code.m
pe, me = PlotkinEncoder(), MatrixEncoder(code)
rng = MersenneTwister(10_010)

pipes = [
    ("RPA",                      me, RPADecoder()),
    ("Dumer min-sum",            pe, DumerDecoder()),
    ("Dumer FHT leaves",         pe, DumerDecoder(leaves = :fht)),
    ("GMD/Dumer-FHT",            pe, GMDDecoder(code, DumerDecoder(leaves = :fht))),
    ("AED-8/Dumer-FHT",          pe, AutomorphismEnsembleDecoder(code, DumerDecoder(leaves = :fht); size = 8, rng)),
    ("DS16 FHT leaves",          pe, DumerShabunovDecoder(16, leaves = :fht)),
    ("GLP cyclic",               pe, GLPDecoder(code, 4 * m; perms = :cyclic)),
    ("GLP pairs (PmCr)",         pe, GLPDecoder(code, 2 * binomial(m, 2); perms = :pairs)),
]

println("Threads: $(nthreads()), trials/point: $TRIALS, Eb/N0 points: $(collect(EBN0_DB))")
println("Resuming RM(2, 10): n=$(blocklength(code)), k=$(dimension(code)), d=$(minimum_distance(code))")
flush(stdout)

for (name, enc, dec) in pipes
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

println("\nRESUME_DONE")
flush(stdout)
