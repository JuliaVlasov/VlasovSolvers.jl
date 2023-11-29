# Quickstart - Landan damping simulation

```@example 1
using VlasovSolvers, Plots

dev = CPU()
```

For now, only the `CPU` device is avalaible but we hope it will be possible 
to run simulation also on `GPU`. 

## Vlasov-Poisson 

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


## Landau damping

- [Landau damping - Wikipedia](https://en.wikipedia.org/wiki/Landau_damping)

We initialize the numerical domain by defining two grids, one for space, one for velocity:

```@docs
VlasovSolvers.OneDGrid
```

```@example 1
nx, nv = 64, 64              # grid resolution
xmin, xmax = 0, 4π           # X Domain length (m)
vmin, vmax = -6, 6           # V Domain length (m)
xgrid = OneDGrid(dev, nx, xmin, xmax)
vgrid = OneDGrid(dev, nv, vmin, vmax)
```

We instantiate a distribution function type to store its values and the grids.

```@docs
VlasovSolvers.DistributionFunction
```

```@example 1
f = DistributionFunction( xgrid, vgrid )
```

We initialize the distribution function for the Landau damping test case:

```@docs
VlasovSolvers.landau!
```

```@example 1
α  = 0.001                   # Perturbation amplitude
kx = 0.5                     # Wave number of perturbation
landau!(f, α, kx)
```

To run the simulation we need to create a [`VlasovProblem`](@ref) and choose a 
method to compute advections. We use the backward semi-lagrangian method with 
5th order b-splines for interpolation [`BSLSpline`](@ref).

```@docs
VlasovSolvers.VlasovProblem
```

```@docs
VlasovSolvers.BSLSpline
```

```@example 1
prob = VlasovProblem(f, BSLSpline(5), dev)
```

```@example 1

stepper = StrangSplitting()  # timestepper
dt = 0.1                     # timestep
nsteps = 1000                # total number of time-steps

```


```@example 1

sol = solve(prob, stepper, dt, nsteps )

t = sol.times

plot(sol; label = "E")
plot!(t, -0.1533*t.-5.50; label="-0.1533t.-5.5")
```

## Simulation with Lagrange interpolation

```@example 1
landau!(f, α, kx)

prob = VlasovProblem(f, BSLLagrange(9), dev)

sol = solve(prob, stepper, dt, nsteps )

t = sol.times

plot(sol; label = "E")
plot!(t, -0.1533*t.-5.50; label="-0.1533t.-5.5")
```
