# Two-stream instability

```@setup tsi
using Plots, Statistics, FFTW, LinearAlgebra
using VlasovSolvers
```

## Parameters
```math
\epsilon = 0.01, \xi = 0.90, v_0 = 2.4
```

## Distribution function

```math
f(x, v) = (1 + \epsilon (( \cos (4 \pi x) + \cos (6 \pi x)) / 1.2 
+ \cos (2 \pi x)))  \frac{1}{\sqrt{2 \pi}}  \frac{2-2 \xi}{3-2 \xi}
(1 + \frac{5 * v^2}{1-\xi})  \exp (- \frac{v^2}{2})
```

## Simulation

```@example tsi
dev = CPU()                  # device
nx, nv = 320, 64             # grid resolution
stepper = StrangSplitting()  # timestepper
dt = 0.1                     # timestep
nsteps = 1000                # total number of time-steps

xmin, xmax = 0, 20π          # X Domain length (m)
vmin, vmax = -8, 8           # V Domain length (m)

xgrid = OneDGrid(dev, nx, xmin, xmax)
vgrid = OneDGrid(dev, nv, vmin, vmax)

df = DistributionFunction( xgrid, vgrid )

two_stream_instability!(df)

contour(vgrid.points, xgrid.points, df.values)
```

```@example tsi
"""
    compute_e(f)

compute Ex using that -ik*Ex = rho 
"""
function compute_e( f )

   dv = f.vgrid.step
   rho = dv * sum(real(f.values), dims=2)
   rho = vec(rho .- mean(rho))
   nx = f.xgrid.len
   xmin = f.xgrid.start
   xmax = f.xgrid.stop
   kx =  2π / (xmax - xmin)
   modes = zeros(Float64, nx)
   modes .= kx * vcat(0:div(nx,2)-1,-div(nx,2):-1)
   modes[1] = 1.0
   rhok = fft(rho) ./ modes
   rhok .*= -1im
   ifft!(rhok)
   real(rhok)

end

import VlasovSolvers: advection!

f = copy(df.values)
fᵗ = transpose(f) |> collect

ex = compute_e(df)

advection!(fᵗ, vgrid, ex, 0.5dt)

dt = 0.01     # Time step
nt = 5000
v = collect(vgrid.points)

anim = @animate for it in 1:nt

    advection!(f, xgrid, v, dt)
    df.values .= f
    ex = compute_e( df )
    transpose!(fᵗ, f)
    advection!(fᵗ, vgrid, ex, dt)
    transpose!(f, fᵗ)
    contourf(vgrid.points, xgrid.points, f, clims=(-1,5))

end every 100

gif(anim, "assets/tsi.gif", fps = 15)
```