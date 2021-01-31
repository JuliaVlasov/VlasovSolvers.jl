# 1D Vlasov–Ampere system

```math
\frac{\partial f}{\partial t} + \upsilon \frac{\partial f}{\partial x}
- E(t,x) \frac{\partial f}{\partial \upsilon} = 0
```

```math
\frac{\partial E}{\partial t} = - J = \int f\upsilon \; d\upsilon
```

## Algorithm 

- For each $j$ compute discrete Fourier transform in $x$ of $f^n(x_i,\upsilon_j)$ yielding $f_k^n(\upsilon_j)$, 

- For $ k \neq 0 $

  - Compute 
     
  ``f^{n+1}_k(\upsilon_j) = e^{−2i\pi k \upsilon \Delta t/L} f_n^k(\upsilon_j),``
     
  - Compute 
     
  ``\rho_k^{n+1} = \Delta \upsilon \sum_j􏰄 f^{n+1}_k(\upsilon_j),``
     
  - Compute
    
  ``E^{n+1}_k = \rho^{n+1}_k L/(2i\pi k \epsilon_0),``
     
  - For $k = 0$ do nothing: 

  ``f_{n+1}(\upsilon_j) = f^n_k(\upsilon_j), E^{n+1}_k = E^n_k``.

  - Perform inverse discrete Fourier transform of $E^{n+1}_k$ and for each $j$ of $f^{n+1}_k (\upsilon_j)$.

```@example 2
using Plots
using VlasovSolvers

dev = CPU()                  # device
nx, nv = 256, 256            # grid resolution
stepper = StrangSplitting()  # timestepper
dt = 0.01                    # timestep
nsteps = 10                  # total number of time-steps

xmin, xmax = 0, 4π           # X Domain length (m)
vmin, vmax = -6, 6           # V Domain length (m)
α  = 0.001                   # Perturbation amplitude
kx = 0.5                     # Wave number of perturbation

xgrid = OneDGrid(dev, nx, xmin, xmax)
vgrid = OneDGrid(dev, nv, vmin, vmax)

f = DistributionFunction( xgrid, vgrid )

landau!(f, α, kx)

prob = VlasovProblem(f, Fourier(xgrid, vgrid), dev)

nsteps = 600
dt = 0.1

sol = solve!(prob, stepper, dt, nsteps )

t =  range(0,stop=60,length=nsteps)

plot(t, -0.1533*t.-5.48)
plot!(t, sol, label="ampere" )
```
