# QuadraticOutputSystems

Collection of tools for working with quadratic output systems defined as
```math
\begin{align*}
\dot x(t) &= Ax(t) + Bx(t), \\
y(t) &= x(t)^{\mathsf{T}}M x(t),
\end{align*}
```
where ``A \in \R^{n \times n}``, ``B \in \R^{n \times m}``, and ``M \in \R^{n \times n}``. As in [BenPD2022](#References) we assume ``M = M^\mathsf{T}``. This can be imposed without loss of generality.

## Functionality

```@docs
h2norm
h2error
h2error_sqr
qo_observability_gramian
```

## References

```latex
@Article{BenPD2022,
  author={Benner, Peter and Goyal, Pawan and Duff, Igor Pontes},
  journal={IEEE Transactions on Automatic Control}, 
  title={Gramians, Energy Functionals, and Balanced Truncation for Linear Dynamical Systems With Quadratic Outputs}, 
  year={2022},
  volume={67},
  number={2},
  pages={886-893},
  doi={10.1109/TAC.2021.3086319}
}
```
