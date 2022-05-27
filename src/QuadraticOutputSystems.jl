module QuadraticOutputSystems

using MatrixEquations, LinearAlgebra

function h2norm(A, B, M)
  Q = qo_observability_gramian(A, B, M)
  return sqrt(tr(B'*Q*B))
end

"""
qo_observability_gramian(A, B, M)

Computes the quadratic-output observability gramian as defined in [1] for a qo system with system matrices ``A, B, M``.
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

export qo_observability_gramian, h2norm

end
