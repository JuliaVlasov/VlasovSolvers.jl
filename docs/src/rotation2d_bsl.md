# Rotation of a gaussian distribution

```math
\frac{df}{dt} +  (y \frac{df}{dx} - x \frac{df}{dy}) = 0
```

```@example rotation2d_bsl
using Plots
using VlasovSolvers
```

```@example rotation2d_bsl

dev = CPU()
n1, n2 = 32, 64
mesh1 = OneDGrid(-π, π, n1)
mesh2 = OneDGrid(-π, π, n2)
x = mesh1.points
y = mesh2.points

tf = 200 * pi
dt = tf/nt

f = DistributionFunction( mesh1, mesh2 )

for (i, xp) in enumerate(x), (j, yp) in enumerate(y)
    xn = cos(tf)*xp - sin(tf)*yp
    yn = sin(tf)*xp + cos(tf)*yp
    f.values[i,j] = exp(-(xn-1)*(xn-1)/0.2)*exp(-(yn-1)*(yn-1)/0.2)
end

```
