using Random, LinearAlgebra

n, m = 10, 3
A = rand(MersenneTwister(0), n, n)
B = rand(MersenneTwister(0), n, m)
M = rand(MersenneTwister(0), n, n)
M = M+M' # M is assumed symmetric.
A = A-I*norm(A)*1.1 # make A stable.

P = QuadraticOutputSystems.controllability_gramian(A, B)
@test norm(A*P+P*A'+B*B') < 1e-8
@test P == P'
@test all(eigvals(P) .>= 0)

Q = qo_observability_gramian(A, B, M)
@test Q == Q'

α = h2norm(A, B, M)
@test α >= 0

