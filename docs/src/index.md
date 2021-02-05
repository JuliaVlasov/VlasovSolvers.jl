# VlasovSolvers.jl Documentation

*First draft of a Vlasov solvers Julia suite.*

## Numerical methods

- Backward Semi-Lagrangian
- Pseudo-Spectral

## Equations

The Vlasov equation can be split into two parts 
```math
\partial_t f + v\cdot \nabla_x f=0 \;\; \mathrm{ and } \;\; \partial_t f +
E\cdot \nabla_v f=0, 
```
so that the numerical solution can be approximated by the successive
solution of each part of this splitting. Indeed, if we denote
``\varphi^x_t(f_0)`` the solution at time $t$ of
```math
\partial_t f + v\cdot \nabla_x f=0, 
```
and $\varphi^v_t(f_0)$ the solution at time $t$ of
```math
\partial_t f + E\cdot \nabla_v f=0, 
```
then the solution $\varphi_t(f_0)$ of the Vlasov equation can
be approximated by
```math
\varphi_t(f_0) \approx \varphi^v_t \circ \varphi^x_t(f_0).  
```
This is the Lie-Trotter splitting which is a first order approximatoin
of the solution $f(t)$. High  order splittings can be derived by
choosing a sequence of coefficients $(a_i, b_i)$ so that
```math
\varphi_t(f_0) \approx \Pi_{i} \;\;  \Big(\varphi^v_{a_i t} \circ \varphi^x_{b_i t}\Big)(f_0).  
```
The splitting enables to reduce the original problem into several
smaller problems which are easier to solve. Indeed, in the Vlasov
case, we are lead to solve one dimensional linear transport equations
which we choose to solve using the semi-Lagrangian method. 
The basics of the semi-Lagrangian method are recalled in the
following. 

We  consider the following linear transport equation 
```math
\partial_t f + a \partial_x f=0, f(t=0, x) = f_0(x), x\in [0, L], 
```
with periodic boundary conditions. We know that the exact solution can
be written as
```math
f(t, x) = f_0(x-at), 
```
or if we consider the solution from  $s$ to $t$`, we
have
```math
f(t, x) = f(s, x-a(t-s)).  
```
The property that $f$ is constant along the characteristics
paves the way of the semi-Lagrangian method. Indeed, let consider
a time discretization $t_n = n\Delta t, n\in \mathbb{N}$,
$\Delta t>0$ and a space discretization
$x_i=i\Delta x, i\in \mathbb{N}, \Delta x>0, \Delta x=L/N_x$
where $N_x$ is the number of points. 
Considering $t=t_{n+1}, s=t_n$ and $x=x_i$ in the former equation, we get 
```math
f(t_{n+1}, x_i) = f(t_n, x_i-a \Delta t).  
```
We assume that $f(t_n, x_i)$ are all known, we have to
compute $f(t_n, x_i-a \Delta t)$ and this is done by a
standard interpolation method. For Vlasov problem, high order
interpolation are required (citer Francis, Sonnen, Michel...). 
Within the framework of this package, there are two possibilities:
Lagrange or B-splines interpolation (of odd degree). 
Arbitrary high order is available and have been tested according to
recent estimates (citer Michel-Latu-Sonnen). For $p$ order Lagrange
interpolation, we have
```math
f(t_n, x) = \sum_{i=0}^k f(t_n, x_i) L_{i,p}(x), 
```
where 
```math
L_{i,p}(x) = \Pi_{0\leq k\leq p, k\neq i} \frac{x-x_k}{x_i-x_k}, \;\;
\mathrm{ for } 0\leq i\leq p.  
```
For the B-spline interpolation of order $p$, we have
```math
f(t_n, x) = \sum_{i=0}^{N_x-1} \eta_i B_{i,p}(x) 
```
where 
```math
B_{i,0}(x) := \left\{
\begin{matrix}
1 & \mathrm{if}  \quad t_i \leq x < t_{i+1} \\
0 & \mathrm{otherwise} 
\end{matrix}
\right.
```
```math
B_{i,p}(x) := \frac{x - t_i}{t_{i+p} - t_i} B_{i,p-1}(x) 
+ \frac{t_{i+p+1} - x}{t_{i+p+1} - t_{i+1}} B_{i+1,p-1}(x).
```
and the coefficients $(\eta_i)_{i=0, \dots, N_x-1}$ are
solution of a linear system to solve.

## Algorithm

Once the linear advection (or transport) equation are solved using the
semi-Lagrangian method, the algorithm to solve the Vlasov-Poisson
equation can be written. To do so, we consider a mesh in the truncated
velocity domain $[-v_{\max}, v_{\max}]$: $v_j = -v_{\max} + j\Delta v,
j=0, \dots, N_v-1, \Delta v=2v_{\max} / N_v$, $N_v$ being the number
of points in the velocity direction. The algorithm based on
Lie-Trotter and semi-Lagrangian method is :

1. Initialization. From the given initial condition $f_0(x, v)$ we can compute the initial electric field $E_0(x)$. 
2. From $t_n$ to $t_{n+1}$. Knowing all the grid point values of $f^n$ and $E^n$
  + Compute $f^\star$ solving
    ```math
     \partial_t f + v\partial_x f=0, 
    ```
using the semi-Lagrangian method $f^\star(x_i, v_j) =
f^n(x_i-v_j\Delta t, v_j)$.
  + Solve the electric field $E^\star$  from the Poisson
equation
    ```math
    \partial_x E^\star = \sum_{j=0}^{N_v-1} f^\star(x_i, v_j) \Delta v -1.  
    ```
  + Compute $f^{n+1}$ solving
    ```math
    \partial_t f + E^\star\partial_v f = 0, 
    ```
    using the semi-Lagrangian method $f^{n+1}(x_i, v_j) = f^n(x_i, v_j-E^\star(x_i) \Delta t)$.
  
## Numerical method for Poisson equation

The Poisson equation is solved using Fourier techniques.
First, we consider the following DFT (Discrete Fourier Transform) of a
$L-$periodic function $g$ defined on a mesh of $N$ points such that
$x_j=j \Delta x, \Delta x=L/N, 0\leq j\leq N-1$ 
```math
\hat{g}_k := \frac{1}{N}\sum_{j=0}^{N-1} g(x_j) e^{-i \frac{2\pi}{L} k
x_j}, \;\; k=0, \dots, N-1. 
```
Thus, from the Poisson equation,  
```math
\partial_x E = \rho -1, 
```
we consider the Fourier transform to get
```math
i k \hat{E}_k = \hat{\rho}_k \;\; \mathrm{  if  } \;\;  k\neq 0, \;\;
\hat{E}_0=0.  
```
Then, we can compute the Fourier coefficients $\hat{E}_k$ as
```math
\hat{E}_k = \frac{1}{ik}\hat{\rho}_k,  
```
and an inverse Fourier transform is used to get $E_i\approx E(x_i)$. 

The extension of this algorithm to the multidimensional case requires
to introduce the electric potential $\varphi$ such that $E=-\nabla_x
\phi$.
The scalar electric potential solves an elliptic Poisson equation
```math
-\Delta \phi = \rho -1
```
which is also solved using Fourier techniques. Now we get (now $k$ is
a vectorial wavenumber) 
```math
|k|^2 \hat{\phi}_k = \hat{\rho}_k \;\; \mathrm{  if  } \;\;  k\neq 0, \;\;
\hat{\phi}_0=0.  
```
so that $\hat{\phi}_k = \hat{\rho}_k/|k|^2$ and the electric field can
be computed using
```math
\hat{E}_k = -ik \hat{\phi}_k  \;\; \mathrm{  if  } \;\;  k\neq 0, \;\;
\hat{E}_0=0.  
```
and an inverse Fourier tranform.

Regarding the semi-Lagrangian method applied to the multidimensional
advection with $x\in \mathbb{R}^d$ and $v\in \mathbb{R}^d$ 
```math
\partial_t f + v\cdot \nabla_x f = 0 \;\; \mathrm{ or } \;\;
\partial_t f + E\cdot \nabla_v f = 0, 
```
we can use the fact that these advection can be split (exactly) into
one-dimensional linear transport equations. Indeed for the
$x$-transport part we have to successively solve 
```math
\partial_t f + v_\alpha\partial_{x_\alpha} f= 0, \;\; \mathrm{ for }\;\;  \alpha=1, \dots, d
```
whereas for the
$v$-transport part we have to successively solve 
```math
\partial_t f + E_\alpha\partial_{v_\alpha} f= 0, \;\; \mathrm{ for } \;\;  \alpha=1, \dots,
d. 
```

## Types and functions

```@autodocs
Modules = [VlasovSolvers]
Order   = [:type, :function]
```
