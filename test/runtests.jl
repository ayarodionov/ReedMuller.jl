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
    @test_throws ArgumentError DumerDecoder(combine = :magic)
    @test_throws ArgumentError BSC(0.7)
    @test_throws ArgumentError BIAWGN(0.0)
end

end
