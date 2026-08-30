using ReedMuller
using Test
using Random

@testset "ReedMuller" begin

@testset "code parameters" begin
    c = RMCode(2, 5)
    @test blocklength(c) == 32
    @test dimension(c) == 16
    @test minimum_distance(c) == 8
    @test rate(c) == 0.5
    @test_throws ArgumentError RMCode(3, 2)
    @test_throws ArgumentError RMCode(-1, 2)
end

@testset "generator matrices" begin
    for (r, m) in [(0, 3), (1, 3), (2, 4), (2, 5), (3, 3), (4, 4)]
        c = RMCode(r, m)
        for b in (:monomial, :plotkin)
            G = generator_matrix(c; basis = b)
            @test size(G) == (dimension(c), blocklength(c))
        end
        # Both bases span the same code: every plotkin row must be
        # orthogonal to the dual code RM(m-r-1, m).
        if r < m
            Gd = generator_matrix(RMCode(m - r - 1, m))
            Gp = generator_matrix(c; basis = :plotkin)
            for i in axes(Gp, 1), j in axes(Gd, 1)
                @test iseven(count(Gp[i, :] .& Gd[j, :]))
            end
        end
    end
end

@testset "encoders agree" begin
    rng = MersenneTwister(1)
    for (r, m) in [(1, 4), (2, 4), (2, 5), (3, 5), (4, 4)]
        c = RMCode(r, m)
        me = MatrixEncoder(c; basis = :plotkin)
        pe = PlotkinEncoder()
        for _ in 1:20
            msg = bitrand(rng, dimension(c))
            @test encode(me, c, msg) == encode(pe, c, msg)
        end
    end
end

@testset "noiseless roundtrip" begin
    rng = MersenneTwister(2)
    cases = [
        (RMCode(1, 4), MatrixEncoder(RMCode(1, 4)), FHTDecoder()),
        (RMCode(1, 6), MatrixEncoder(RMCode(1, 6)), FHTDecoder()),
        (RMCode(2, 5), MatrixEncoder(RMCode(2, 5)), ReedDecoder()),
        (RMCode(0, 4), MatrixEncoder(RMCode(0, 4)), ReedDecoder()),
        (RMCode(4, 4), MatrixEncoder(RMCode(4, 4)), ReedDecoder()),
        (RMCode(2, 5), PlotkinEncoder(), DumerDecoder()),
        (RMCode(3, 6), PlotkinEncoder(), DumerDecoder(combine = :exact)),
        (RMCode(4, 4), PlotkinEncoder(), DumerDecoder()),
        (RMCode(0, 5), PlotkinEncoder(), DumerDecoder()),
        (RMCode(2, 5), PlotkinEncoder(), DumerShabunovDecoder(4)),
        (RMCode(3, 6), PlotkinEncoder(), DumerShabunovDecoder(8, combine = :exact)),
        (RMCode(4, 4), PlotkinEncoder(), DumerShabunovDecoder(2)),
        (RMCode(0, 5), PlotkinEncoder(), DumerShabunovDecoder(1)),
        (RMCode(2, 4), MatrixEncoder(RMCode(2, 4)), SidelnikovPershakovDecoder()),
        (RMCode(2, 6), MatrixEncoder(RMCode(2, 6)), SidelnikovPershakovDecoder()),
    ]
    for (c, enc, dec) in cases
        for _ in 1:20
            msg = bitrand(rng, dimension(c))
            cw = encode(enc, c, msg)
            @test decode(dec, c, hard_llr(cw)) == msg
        end
    end
end

@testset "Reed corrects up to (d-1)/2 errors" begin
    rng = MersenneTwister(3)
    for (r, m) in [(1, 4), (2, 5)]
        c = RMCode(r, m)
        enc = MatrixEncoder(c)
        t = (minimum_distance(c) - 1) ÷ 2
        for _ in 1:50
            msg = bitrand(rng, dimension(c))
            cw = encode(enc, c, msg)
            pos = randperm(rng, blocklength(c))[1:t]
            y = copy(cw)
            y[pos] .⊻= true
            @test decode(ReedDecoder(), c, hard_llr(y)) == msg
        end
    end
end

@testset "FHT is maximum likelihood" begin
    rng = MersenneTwister(4)
    c = RMCode(1, 4)
    enc = MatrixEncoder(c)
    k, n = dimension(c), blocklength(c)
    codebook = [encode(enc, c, BitVector(digits(Bool, w, base = 2, pad = k)))
                for w in 0:(2^k - 1)]
    ch = BIAWGN(1.0)
    for _ in 1:50
        msg = bitrand(rng, k)
        llr = transmit(rng, ch, encode(enc, c, msg))
        est_cw = encode(enc, c, decode(FHTDecoder(), c, llr))
        # ML codeword maximises correlation Σ llr_i (1 - 2 c_i)
        score(cw) = sum(llr[i] * (cw[i] ? -1.0 : 1.0) for i in 1:n)
        @test score(est_cw) ≈ maximum(score, codebook)
    end
end

@testset "Dumer-Shabunov list decoding" begin
    rng = MersenneTwister(6)
    c = RMCode(2, 5)
    enc = PlotkinEncoder()
    ch = BIAWGN(0.9)

    # With L = 1 the list decoder must agree with plain Dumer decoding.
    for _ in 1:50
        llr = transmit(rng, ch, encode(enc, c, bitrand(rng, dimension(c))))
        @test decode(DumerShabunovDecoder(1), c, llr) == decode(DumerDecoder(), c, llr)
    end

    # A larger list must beat plain Dumer over a noisy batch.
    errs = Dict(1 => 0, 16 => 0)
    for _ in 1:300
        msg = bitrand(rng, dimension(c))
        llr = transmit(rng, ch, encode(enc, c, msg))
        for L in (1, 16)
            errs[L] += decode(DumerShabunovDecoder(L), c, llr) != msg
        end
    end
    @test errs[16] < errs[1]

    # With the list covering the whole codebook, decoding is exact ML
    # by correlation on a small code.
    small = RMCode(1, 3)                # k = 4, 16 codewords
    k, n = dimension(small), blocklength(small)
    codebook = [encode(enc, small, BitVector(digits(Bool, w, base = 2, pad = k)))
                for w in 0:(2^k - 1)]
    dec = DumerShabunovDecoder(2^k)
    for _ in 1:50
        llr = transmit(rng, BIAWGN(1.2), encode(enc, small, bitrand(rng, k)))
        est_cw = encode(enc, small, decode(dec, small, llr))
        score(cw) = sum(llr[i] * (cw[i] ? -1.0 : 1.0) for i in 1:n)
        @test score(est_cw) ≈ maximum(score, codebook)
    end
end

@testset "Sidelnikov-Pershakov decoding" begin
    rng = MersenneTwister(7)

    # Guaranteed radius: any pattern of up to 2^(m-3) - 1 errors.
    for m in (5, 6)
        c = RMCode(2, m)
        enc = MatrixEncoder(c)
        t = (minimum_distance(c) - 1) ÷ 2
        for _ in 1:30
            msg = bitrand(rng, dimension(c))
            y = encode(enc, c, msg)
            y[randperm(rng, blocklength(c))[1:t]] .⊻= true
            @test decode(SidelnikovPershakovDecoder(), c, hard_llr(y)) == msg
        end
    end

    # Beyond half distance: corrects most random patterns of 1.5t errors.
    c = RMCode(2, 7)
    enc = MatrixEncoder(c)
    t15 = 3 * (minimum_distance(c) ÷ 2) ÷ 2
    ok = 0
    for _ in 1:50
        msg = bitrand(rng, dimension(c))
        y = encode(enc, c, msg)
        y[randperm(rng, blocklength(c))[1:t15]] .⊻= true
        ok += decode(SidelnikovPershakovDecoder(), c, hard_llr(y)) == msg
    end
    @test ok >= 45
end

@testset "simulate" begin
    rng = MersenneTwister(5)
    c = RMCode(1, 5)
    res = simulate(MatrixEncoder(c), FHTDecoder(), c, BSC(0.01);
                   trials = 200, rng)
    @test res.wer <= 0.05
    @test res.ber <= res.wer
    @test_throws ArgumentError simulate(PlotkinEncoder(), FHTDecoder(), c, BSC(0.01))
end

@testset "input validation" begin
    c = RMCode(2, 4)
    @test_throws DimensionMismatch encode(PlotkinEncoder(), c, falses(3))
    @test_throws DimensionMismatch decode(DumerDecoder(), c, zeros(5))
    @test_throws ArgumentError decode(FHTDecoder(), c, zeros(16))
    @test_throws ArgumentError decode(SidelnikovPershakovDecoder(), RMCode(1, 4), zeros(16))
    @test_throws ArgumentError DumerDecoder(combine = :magic)
    @test_throws ArgumentError DumerShabunovDecoder(0)
    @test_throws ArgumentError DumerShabunovDecoder(4, combine = :magic)
    @test_throws ArgumentError BSC(0.7)
    @test_throws ArgumentError BIAWGN(0.0)
end

end
