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


```@setup 1
using Plots
```

## Input parameters

```@example 1
using VlasovSolvers

dev = CPU()                  # device
nx, nv = 64, 64              # grid resolution
stepper = StrangSplitting()  # timestepper
dt = 0.1                     # timestep
nsteps = 1000                # total number of time-steps

xmin, xmax = 0, 4π           # X Domain length (m)
vmin, vmax = -6, 6           # V Domain length (m)
α  = 0.001                   # Perturbation amplitude
kx = 0.5                     # Wave number of perturbation

xgrid = OneDGrid(dev, nx, xmin, xmax)
vgrid = OneDGrid(dev, nv, vmin, vmax)

f = DistributionFunction( xgrid, vgrid )

landau!(f, α, kx)

prob = VlasovProblem(f, BSL(5), dev)

sol = solve!(prob, stepper, dt, nsteps )

t = LinRange(0,100,1000)

plot( t, sol; label = "E")
plot!(t, -0.1533*t.-5.50; label="-0.1533t.-5.5")
```

## Types and functions

```@autodocs
Modules = [VlasovSolvers]
Order   = [:type, :function]
```
