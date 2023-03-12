module QuadraticOutputSystems

using LinearAlgebra
using MatrixEquations

sym(A) = 0.5 * (A + A')

"""
    h2norm(A, B, M)

Computes the h2norm of quadratic output system
``
\\dot x = Ax+Bu, y  = x^{\\mathsf{T}} M x
``
"""
function h2norm(A, B, M)
    Q = qo_observability_gramian(A, B, M)
    return sqrt(tr(B' * Q * B))
end

"""
    h2error(A1, B1, M1, A2, B2, M2)

Computes the h2-error between two quadratic output systems defined by `A1`,
`B1`, and `M1` and `A2`, `B2`, and `M2`, respectively.
Note that the optional named argumentd `BQ1B`, which defaults to
`B1'*qo_observability_gramian(A1, B1, M1)*B1` can be passed, if `BQ1B` is known
beforehand.
"""
function h2error(
    A1,
    B1,
    M1,
    A2,
    B2,
    M2;
    BQ1B = B1' * qo_observability_gramian(A1, B1, M1) * B1,
    X = sylvc(A1, A2', -B1 * B2'),
)
    return sqrt(abs(h2error_sqr(A1, B1, M1, A2, B2, M2; BQ1B, X)))
end

"""
    h2error_sqr(A1, B1, M1, A2, B2, M2)

Computes the **square** of the h2-error between two quadratic output systems
defined by `A1`, `B1`, and `M1` and `A2`, `B2`, and `M2`, respectively.
Note that the optional named argumentd `BQ1B`, which defaults to
`B1'*qo_observability_gramian(A1, B1, M1)*B1` can be passed, if `BQ1B` is known
beforehand.
"""
function h2error_sqr(
    A1,
    B1,
    M1,
    A2,
    B2,
    M2;
    BQ1B=B1' * qo_observability_gramian(A1, B1, M1) * B1,
    X=sylvc(A1, A2', -B1 * B2'),
)
    M1XM2 = M1 * X * M2
    Z = sylvc(A1', A2, -M1XM2)
    Q2 = qo_observability_gramian(A2, B2, M2)
    return tr(BQ1B + B2' * Q2 * B2 - 2 * B1' * Z * B2)
end

"""
    h2error_sqr_fast(A1, B1, M1, A2, B2, M2)

Computes the **square** of the h2-error between two quadratic output systems
defined by `A1`, `B1`, and `M1` and `A2`, `B2`, and `M2`, respectively. This is
faster, if only `M2` is changed and size(A1, 1) << size(A2, 1) because all
larger Sylvester and Lyapunov equations are solved beforehand.
"""
function h2error_sqr_fast(
    A1,
    B1,
    M1,
    A2,
    B2,
    M2;
    BQ1B = B1' * qo_observability_gramian(A1, B1, M1) * B1,
    X = sylvc(A1, A2', -B1 * B2'),
    J2 = sylvester_helper(A1, A2, B1, B2, M1, X),
)
    m = size(B1, 2)
    TwiceB1TZB = reshape(J2 * vec(M2), m, m)
    Q2 = qo_observability_gramian(A2, B2, M2)
    return tr(BQ1B + B2' * Q2 * B2 - TwiceB1TZB)
end

function sylvester_helper(A1, A2, B1, B2, M1, X)
    n = size(A1, 1)
    r = size(A2, 1)
    AS = kron(I(r), A1') + kron(A2', I(n))
    J2h = AS\(kron(I(r), M1 * X))
    J2 = -2*kron(B2', B1')*J2h
    return J2
end

"""
    qo_observability_gramian(A, B, M)

Computes the quadratic-output observability gramian as defined in BenPD2022 for a quadratic output system with system matrices ``A, B, M``.
"""
function qo_observability_gramian(A, B, M)
    @assert M == M' "Please provide a symmetric M. Note that this can be assumed for QO-systems without loss of generality"
    PU = controllability_gramian_factor(A, B)
    U = plyapc(Array(A'), M*PU)
    return U*U'
end

"""
    controllability_gramian(A, B)

Computes the controllability gramian of a system ``\\dot x = Ax+Bu``, which is given by the solution ``P`` of the Lyapunov equation: ``AP+PA^{\\mathsf{T}}+BB^{\\mathsf{T}} = 0``.
Note that we use the function `plyapc` to ensure a positive (semi) definite solution `P >= 0`.
"""
function controllability_gramian(A, B)
    P = plyapc(A, B)
    return P * P'
end

function controllability_gramian_factor(A, B)
    U = plyapc(A, B)
    return U
end

export qo_observability_gramian, h2norm, h2error, h2error_sqr


end
