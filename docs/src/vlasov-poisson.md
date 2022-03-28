# Vlasov-Poisson 

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


## Input parameters

```@example 1
using VlasovSolvers, Plots

dev = CPU()                  # device
nx, nv = 64, 64              # grid resolution
stepper = StrangSplitting()  # timestepper
dt = 0.1                     # timestep
nsteps = 1000                # total number of time-steps

xmin, xmax = 0, 4π           # X Domain length (m)
vmin, vmax = -6, 6           # V Domain length (m)
α  = 0.001                   # Perturbation amplitude
kx = 0.5                     # Wave number of perturbation
```

## Simulation

```@example 1
xgrid = OneDGrid(dev, nx, xmin, xmax)
vgrid = OneDGrid(dev, nv, vmin, vmax)

f = DistributionFunction( xgrid, vgrid )

landau!(f, α, kx)

prob = VlasovProblem(f, BSLSpline(5), dev)

sol = solve(prob, stepper, dt, nsteps )

t = sol.times

plot(sol; label = "E")
plot!(t, -0.1533*t.-5.50; label="-0.1533t.-5.5")
```
