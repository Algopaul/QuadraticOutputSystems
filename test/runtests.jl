using QuadraticOutputSystems
using Test
using Random, LinearAlgebra

function generate_stable_qo_system(n, m)
    A = rand(MersenneTwister(0), n, n)
    B = rand(MersenneTwister(0), n, m)
    M = rand(MersenneTwister(0), n, n)
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
        A, B, _ = generate_stable_qo_system(10, 5)
        P = QuadraticOutputSystems.controllability_gramian(A, B)
        @test norm(A * P + P * A' + B * B') < 1e-8
        @test is_positive_definite(P)
    end

    @testset "QO-Gramian" begin
        A, B, M = generate_stable_qo_system(10, 2)
        Q = qo_observability_gramian(A, B, M)
        @test is_positive_definite(Q)
        # TODO: Check that Q is in fact the QO-gramian.
    end

    @testset "QH2-norm" begin
        A, B, M = generate_stable_qo_system(10, 2)
        α = h2norm(A, B, M)
        @test α >= 0
        M .= 0
        α = h2norm(A, B, M)
        @test α ≈ 0
    end

    @testset "QH2-error" begin
        A1, B1, M1 = generate_stable_qo_system(10, 2)
        A2, B2, M2 = generate_stable_qo_system(15, 2)
        β = h2error(A1, B1, M1, A2, B2, M2)
        @test β >= 0
        β = h2error(A1, B1, M1, A1, B1, M1)
        @test β ≈ 0
    end
end
