using Random, LinearAlgebra

function generate_stable_qo_system(n, m)
  A = rand(MersenneTwister(0), n, n)
  B = rand(MersenneTwister(0), n, m)
  M = rand(MersenneTwister(0), n, n)
  M = M+M' # M is assumed symmetric.
  A = A-I*norm(A)*1.1 # make A stable.
  return A, B, M
end

A, B, M = generate_stable_qo_system(10, 2)

P = QuadraticOutputSystems.controllability_gramian(A, B)
@test norm(A*P+P*A'+B*B') < 1e-8
@test norm(P - P') < 1e-8
@test all(eigvals(P) .>= 0)

Q = qo_observability_gramian(A, B, M)
@test Q == Q'

α = h2norm(A, B, M)
@test α >= 0

A1, B1, M1 = generate_stable_qo_system(10, 2)
A2, B2, M2 = generate_stable_qo_system(15, 2)
β = h2error(A1, B1, M1, A2, B2, M2)
@test α >= 0
