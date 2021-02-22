# VlasovSolvers.jl Documentation

*First draft of a Vlasov solvers Julia suite.*


The Vlasov equation is of paramount importance in plasma physics. It is a nonlinear transport equation 
satisfied by the distribution function $f$ of the considered charged particles (ions or electrons typically).  
The Vlasov equation is a kinetic equation which means the unknown $f$ not only depends on time $t$ and space $x$ 
but also on the velocity repartition of the particles through the variable $v$. Hence, the distribution function 
$f$ depends on $t\geq 0$, $x\in \Omega\subset \mathbb{R}^d$ and $v\in \mathbb{R}^d$ with $d\geq 1$ the dimension of the problem. 

The main goal of this library is to propose some efficient numerical tools to solve numerically the Vlasov-Poisson equation 
using the semi-Lagrangian method. 

We consider a population of electrons whereas the ions are supposed stationary with a constant density equal to one. 
The spatial domain $\Omega$ is a torus in dimension $d$ so that $x\in \Omega=[0, L]^d$ with periodic boundary conditions.  
In this framework, the Vlasov-Poisson equation we intend to solve can be written as  
```math
\partial_t f + v\cdot \nabla_x f + E\cdot \nabla_v f =0, 
```
where the electric field $E=E(t, x)\in \mathbb{R}^d$ depends on the solution $f$ through the Poisson equation 
```math
\nabla_x \cdot E=  \int f dv - 1. 
```
The electric field is supposed to derive from an electric potential $\phi(t, x)\in \mathbb{R}$ so that 
$E(t, x)=-\nabla_x \phi(t, x)$ and $\phi$ solves the following elliptic equation 
```math
-\Delta \phi=  \int f dv - 1. 
```
In the following, some elements detailing the numerical method developed in the library are given 
(splitting, semi-Lagrangian method, Poisson equation solver). 


## Splitting

Splitting methods enable to reduce the resolution of complex systems into the successive  
resolution of simple problems. The way complex systems are split strongly depends on the 
considered problem. For the Vlasov-Poisson system, we appeal to its specific Hamiltonian structure (see Casas-Crouseilles-Faou-Mehrenberger), 
which enables to split the Vlasov equation into the two following parts 
```math
\partial_t f + v\cdot \nabla_x f=0 \;\; \mathrm{ and } \;\; \partial_t f +
E\cdot \nabla_v f=0, 
```
so that the numerical solution can be approximated by the successive
solution of each part. Indeed, if we denote
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

A very interesting aspect of this splitting comes from the fact that each subpart $\varphi^x_t(f_0)$ and $\varphi^v_t(f_0)$ 
can be solved exactly in time. Indeed, the solution at time $t$ of $\partial_t f + v\cdot \nabla_x f=0$ 
with a given initial condition $f(0, x, v)=f_0(x, v)$ is $f(t, x, v)=f_0(x-vt, v)$. 
For a given electric field $E(t, x)$ (typically obtained from the solution of the Poisson equation), the  
solution at time $t$ of $\partial_t f + E\cdot \nabla_v f=0$ 
with a given initial condition $f(0, x, v)=f_0(x, v)$ is $f(t, x, v)=f_0(x, v-E(0, x) t)$. Indeed, it is worth mentioning 
that the electric field $E$ does not depend on time during this step. Since $E$ only depends on $\int f dv$, 
the velocity integration of $\partial_t f + E\cdot \nabla_v f=0$ leads to $\frac{d}{dt}\int f dv = 0$ so that 
$\frac{d}{dt} E(t, x) = 0$ and $E(t, x)=E(0, x)$ along this step. 


This is the Lie-Trotter splitting which is a first order approximatoin
of the solution $f(t)$. High order splittings can be derived by
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

We are faced with multidimensional linear transport equation which can be split again 
into one dimensional linear transport equation. Indeed, for the part 
$\partial_t f + v\cdot \nabla_x f = 0$, we can  split exactly this $d$-dimensional transport equation into 
$$
\partial_t + v_\alpha \partial_{x_\alpha} f = 0, \;\; \alpha=1, \dots, d,   
$$
where $x_\alpha$ (resp. $v_\alpha$) denotes the $\alpha$-th component of $x$ (resp. $v$). 
Similarly, we can  split exactly the part $\partial_t f + E\cdot \nabla_v f = 0$ 
into $d$ one dimensional linear transport equations 
$$
\partial_t + E_\alpha \partial_{v_\alpha} f = 0, \;\; \alpha=1, \dots, d,    
$$
where $E_\alpha$ denotes the $\alpha$-th compoent of the electric field $E$. 

## Semi-Lagrangian method 
According to the previous section, we are lead to solve the following linear transport equation 
```math
\partial_t f + a \partial_x f=0, f(t=0, x) = f_0(x), x\in [0, L], 
```
with periodic boundary conditions. According to the above notations, 
$x$ denotes here the spatial or velocity direction. 

We know that the exact solution of $\partial_t f + a \partial_x f=0$ can be written as
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
$x_i=i\Delta x, i\in N_x, \Delta x>0, \Delta x=L/N_x$
where $N_x$ is the number of points. 
Considering $t=t_{n+1}, s=t_n$ and $x=x_i$ in the former equation, we get 
```math
f(t_{n+1}, x_i) = f(t_n, x_i-a \Delta t).  
```
We assume that $f(t_n, x_i)$ are all known, we have to
compute $f(t_n, x_i-a \Delta t)$ and this is done by a
standard interpolation method. For Vlasov problem, high order
interpolation are required (citer Francis, Sonnen, Michel...). 
Within the framework of this package, we have chosen 
piecewise polynomial interpolation of two kinds:
Lagrange interpolation or B-splines interpolation (of odd degree). 
Arbitrary high order of Lagrange and B-splines polynomials are 
available and have been tested according to recent estimates (citer Michel-Latu-Sonnen). 

For $p$ order Lagrange interpolation, we have
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
solution of a linear system to solve (citer De Boor, Michel, ...).


## Numerical method for Poisson equation

The Poisson equation is solved using Fourier techniques.
First, we consider the following DFT (Discrete Fourier Transform) of a
$L-$periodic function $g$ defined on a mesh of $N_x$ points such that
$x_j=j \Delta x, \Delta x=L/N_x, 0\leq j\leq N_x-1$ 
```math
\hat{g}_k := \frac{1}{N_x}\sum_{j=0}^{N_x-1} g(x_j) e^{-i \frac{2\pi}{L} k
x_j}, \;\; k=0, \dots, N_x-1. 
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
and an inverse Fourier tranform. See book Sonnen for more details. 

## Algorithm

Once the linear advection (or transport) equation are solved using the
semi-Lagrangian method, the algorithm to solve the Vlasov-Poisson
equation can now be written. To do so, we consider a uniform mesh of the phase space $(x, v)$. 
For the spatial domain $[0, L]$, we introduce $x_i=i\Delta x, i=0, \dots, N_x-1$, 
$N_x$ being the number of points in $x$ and $\Delta x=L/N_x$. 
For the velocity direction, we need to consider a truncated velocity domain $[-v_{\max}, v_{\max}]$ 
and the mesh is defined by $v_j = -v_{\max} + j\Delta v,
j=0, \dots, N_v-1, \Delta v=2v_{\max} / N_v$, $N_v$ being the number
of points in the velocity direction. 

Then, the algorithm based on Lie-Trotter and semi-Lagrangian method is:

1. Initialization. From the given initial condition $f_0$, we evaluate it on the phase space mesh to get $f^0_{i,j}=f_0(x_i, v_j)$. We then compute the initial electric field $E_{0, i}\approx E_0(x_i)$. 
2. From $t_n$ to $t_{n+1}$. Knowing all the grid point values of $f^n$ and $E^n$
  + Compute $f^\star_{i,j}$ solving
```math
\partial_t f + v\partial_x f=0, 
```
using the semi-Lagrangian method $f^\star_{i, j} \approx f^n(x_i-v_j\Delta t, v_j)$.
  + Solve the electric field $E^\star_i$  from the Poisson
equation
```math
\partial_x E^\star = \sum_{j=0}^{N_v-1} f^\star_{i, j} \Delta v -1.  
```
  + Compute $f^{n+1}$ solving
```math
\partial_t f + E^\star\partial_v f = 0, 
```
using the semi-Lagrangian method $f^{n+1}_{i, j} approx f^n(x_i, v_j-E^\star_i \Delta t)$. 

The extension to the well-known Strang  splitting (which is second order accurate in time) can be written as follows. 

1. Initialization. From the given initial condition $f_0$, we evaluate it on the phase space mesh to get $f^0_{i,j}=f_0(x_i, v_j)$. We then compute the initial electric field $E_{0, i}\approx E_0(x_i)$. 
  
2. From $t_n$ to $t_{n+1}$. Knowing all the grid point values of $f^n_{i,j}$ and $E^n_i$
  + Compute $f^\star_{i,j}$ solving
```math
\partial_t f + v\partial_x f=0, 
```
using the semi-Lagrangian method $f^\star_{i, j} \approx f^n(x_i-v_j\Delta t/2, v_j)$.
  + Solve the electric field $E^\star_i$  from the Poisson
equation
```math
\partial_x E^\star = \sum_{j=0}^{N_v-1} f^\star_{i, j} \Delta v -1.  
```
  + Compute $f^{\star\star}$ solving
```math
\partial_t f + E^\star\partial_v f = 0, 
```
using the semi-Lagrangian method $f^{\star\star}_{i, j} approx f^n(x_i, v_j-E^\star_i \Delta t)$. 
  
+ Compute $f^\star_{i,j}$ solving
```math
\partial_t f + v\partial_x f=0, 
```
using the semi-Lagrangian method $f^{n+1}{i, j} \approx f^{\star\star}(x_i-v_j\Delta t/2, v_j)$.

