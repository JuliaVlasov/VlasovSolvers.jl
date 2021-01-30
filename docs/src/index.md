# VlasovSolvers.jl Documentation

First draft of a Vlasov solver.

## Vlasov-Poisson equation

We consider the dimensionless Vlasov-Poisson equation for one species
with a neutralizing background.

```math
\frac{\partial f}{\partial t}+ v\cdot \nabla_x f + E(t,x) \cdot \nabla_v f = 0,
```

```math
- \Delta \phi = 1 - \rho, E = - \nabla \phi
```

```math
\rho(t,x)  =  \int f(t,x,v)dv.
```

- [Vlasov Equation - Wikipedia](https://en.wikipedia.org/wiki/Vlasov_equation)


```@autodocs
Modules = [VlasovSolvers]
Order   = [:type, :function]
```
