using QuadraticOutputSystems
using Test
using Random, LinearAlgebra

function generate_stable_qo_system!(rng, n, m)
    A = rand(rng, n, n)
    B = rand(rng, n, m)
    M = rand(rng, n, n)
    M = M + M' # M is assumed symmetric.
    A = A - I * norm(A) * 1.1 # make A stable.
    return A, B, M
end

function is_positive_definite(A)
    try
        cholesky(A)
        return true
    catch
        return false
    end
end


@testset "QuadraticOutputSystems.jl" verbose = true begin

    @testset "Controllability Gramian" begin
        A, B, _ = generate_stable_qo_system!(MersenneTwister(0), 10, 5)
        P = QuadraticOutputSystems.controllability_gramian(A, B)
        @test norm(A * P + P * A' + B * B') < 1e-8
        @test is_positive_definite(P)
    end

    @testset "QO-Gramian" begin
        A, B, M = generate_stable_qo_system!(MersenneTwister(0), 10, 2)
        Q = qo_observability_gramian(A, B, M)
        @test is_positive_definite(Q)
        P = QuadraticOutputSystems.controllability_gramian(A, B)
        @test norm(A'*Q + Q*A + M*P*M) < 1e-6
    end

    @testset "QH2-norm" begin
        A, B, M = generate_stable_qo_system!(MersenneTwister(0), 10, 2)
        α = h2norm(A, B, M)
        @test α >= 0
        M .= 0
        α = h2norm(A, B, M)
        @test α ≈ 0
    end

    @testset "QH2-error" begin
        rng = MersenneTwister(0)
        A1, B1, M1 = generate_stable_qo_system!(rng, 25, 2)
        A2, B2, M2 = generate_stable_qo_system!(rng, 15, 2)
        β = h2error(A1, B1, M1, A2, B2, M2)
        @test β >= 0
        β = h2error(A1, B1, M1, A1, B1, M1)
        @test β < 1e-5
    end

    @testset "QH2-error-fast" begin
        rng = MersenneTwister(0)
        A1, B1, M1 = generate_stable_qo_system!(rng, 20, 2)
        A2, B2, M2 = generate_stable_qo_system!(rng, 15, 2)
        β1 = h2error_sqr(A1, B1, M1, A2, B2, M2)
        β2 = QuadraticOutputSystems.h2error_sqr_fast(A1, B1, M1, A2, B2, M2)
        @test β1 ≈ β2
        β = QuadraticOutputSystems.h2error_sqr_fast(A1, B1, M1, A1, B1, M1)
        @test β < 1e-5
    end
end
