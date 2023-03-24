# gpu_finite_element
GPU accelerated structural stress analysis using the finite element method
or: [building bridges](https://www.youtube.com/watch?v=tVLxcZMzZmw) in the
computer for fun.

Compile with `make gpu`. Has non-GPU version too without CUDA dependency (`make` or `make openmp`).

## Notes
If $\psi^j$ is a small displacement field with $j$ indexing the spatial directions, then the strain is
$$
s^{ij} = \frac{1}{2}(\nabla^i \psi^j + \nabla^j \psi^i)
$$
and for isotropic linear elasticity (scalar $E$), the stress is
$$
\sigma^{ij} = E s^{ij}
$$
and Newton's second law $(F = ma)$ with density $\rho$ and external forces per volume $f^j$ is
$$
\rho \frac{\partial^2 \psi^j}{\partial t^2} = \nabla_i \sigma^{ij} + f^j
$$

To solve for the steady state ($\partial_t^2 \psi^j = 0$), discretize by choosing some
basis functions $u(x)$ (representing the degrees of freedom of interest) and expanding
$$
\psi^j(x) = \sum_n c_n u_n^j(x)
$$
and then require that the coefficients $c$ make the following inner product zero with each $u$ (weak form discretized):
$$
0 = \int u_{mj} (\nabla_i \sigma^{ij} + f^j)\,d^3x
$$

This gives the matrix problem
$$
\sum_{n} A_{mn} c_n = b_m
$$
where
$$
A_{mn} = -\frac{E}{2}\int u_{mj} \nabla_i (\nabla^i u_n^j + \nabla^j u_n^i)\,d^3x\\
b_m = \int u_{mj}f^j\,d^3x
$$
and $A$ (the stiffness matrix) is symmetric (and sparse if $u_m$ was chosen well).

The conjugate gradient method can iteratively solve this for the coefficients
$c$ and give the discretized solution.
