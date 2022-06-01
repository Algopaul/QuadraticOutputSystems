module QuadraticOutputSystems

using MatrixEquations, LinearAlgebra

"""
`h2norm(A, B, M)`

Computes the h2norm of quadratic output system
``
\\dot x = Ax+Bu, y  = x^{\\mathsf{T}} M x
``
"""
function h2norm(A, B, M)
  Q = qo_observability_gramian(A, B, M)
  return sqrt(tr(B'*Q*B))
end

"""
`h2error(A1, B1, M1, A2, B2, M2)`

Computes the h2-error between two quadratic output systems defined by `A1`,
`B1`, and `M1` and `A2`, `B2`, and `M2`, respectively.
Note that the optional named argumentd `BQ1B`, which defaults to
`B1'*qo_observability_gramian(A1, B1, M1)*B1` can be passed, if `BQ1B` is known
beforehand.
"""
function h2error(A1, B1, M1, A2, B2, M2; BQ1B = B1'*qo_observability_gramian(A1, B1, M1)*B1)
  X = sylvc(A1, A2', -B1*B2')
  M1XM2 = M1*X*M2
  Z = sylvc(A1', A2, -M1XM2)
  Q2 = qo_observability_gramian(A2, B2, M2)
  return sqrt(tr(BQ1B+B2'*Q2*B2-2*B1'*Z*B2))
end

"""
`qo_observability_gramian(A, B, M)`

Computes the quadratic-output observability gramian as defined in BenPD2022 for a quadratic output system with system matrices ``A, B, M``.
"""
function qo_observability_gramian(A, B, M)
  @assert M == M' "Please provide a symmetric M. Note that this can be assumed for QO-systems without loss of generality"
  P = controllability_gramian(A, B)
  MPM = M*P*M
  MPM = 0.5*(MPM+MPM')
  Q = lyapc(Array(A'), M*P*M)
  Q = 0.5*(Q+Q')
end

"""
controllability_gramian(A, B)

Computes the controllability gramian of a system ``\\dot x = Ax+Bu``, which is given by the solution ``P`` of the Lyapunov equation: ``AP+PA^{\\mathsf{T}}+BB^{\\mathsf{T}} = 0``.
Note that we use the function `plyapc` to ensure a positive (semi) definite solution `P >= 0`.
"""
function controllability_gramian(A, B)
  P = plyapc(A, B)
  return P*P'
end

export qo_observability_gramian, h2norm, h2error

end
