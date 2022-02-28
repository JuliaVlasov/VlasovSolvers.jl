# Vlasov HMF

In this example we will show how the VlasovSolvers package can be used to solve the Vlasov equation in the Hamiltonian Mean Field framework (usually called Vlasov-HMF).

The Vlasov-HMF system is a simplification of the Vlasov-Poisson system, which is itself a simplification of the Vlasov-Maxwell equations.


Recall the Vlasov-Poisson equations in the noncollisional case:
```math
\partial_t f(t,x,v) + v\cdot \partial_x f(t,x,v) + E(t,x)\cdot \partial_v f(t,x,v) = 0
```

The quantity $E(t,x)$ is the electric field, defined in the Poisson framework by
```math
-\Delta \Phi = 1 - \rho,\, E = -\nabla \Phi,\, \rho(t,x) = \int f(t,x,v)dv
```

The Poisson equation $-\Delta \Phi = 1 - \rho$ on a periodic space-domain $[0, L]$ is usually solved by means of a Fourier
transform, assuming $\Phi$ has a zero average. In the discrete case this is performed by a DFT, involving $N_x$ Fourier
 modes in total (where $N_x$ is the number of spatial nodes). In the HMF framework, we apply the same idea but restrict
ourselves to the Fourier modes corresponding to the frequencies $k=\pm 1$ (the mode corresponding to $k=0$ is zero 
since $\Phi$ has a zero average).

```@example hmf
using LinearAlgebra, QuadGK, Roots, FFTW
using VlasovSolvers
using Plots
```

```@example hmf
function mag(β, mass)
    F(m) = begin
        g(x, n, m) = (1 / π) * (exp(β * m * cos(x)) * cos(n * x))
        bessel0(x) = g(x, 0, m)
        bessel1(x) = g(x, 1, m)
        mass * quadgk(bessel1, 0, π)[1] / quadgk(bessel0, 0, π)[1] - m
    end
    find_zero(F, (0, mass))
end
```

```@example hmf
function Norm(f::Array{Float64,2}, delta1, delta2)
    delta1 * sum(delta2 * sum(real(f), dims=1))
end
```

```@example hmf
"""
    hmf_poisson!(fᵗ    :: Array{Complex{Float64},2},
                 mesh1 :: OneDGrid,
                 mesh2 :: OneDGrid,
                 ex    :: Array{Float64})

    Compute the electric hamiltonian mean field from the
    transposed distribution function

"""
function hmf_poisson!(fᵗ::Array{Complex{Float64},2},
        mesh1::OneDGrid,
        mesh2::OneDGrid,
        ex::Array{Float64})

    n1 = mesh1.len
    rho = mesh2.step .* vec(sum(fᵗ, dims=1)) # ≈ ∫ f(t,x_i,v)dv, i=1, ..., N_x
    kernel = zeros(Float64, n1)
    k = -(mesh1.stop - mesh1.start) / (2π)
    kernel[2]   =  k    # fourier mode  1
    kernel[end] = -k    # fourier mode -1
    ex .= real(ifft(fft(rho) .* 1im .* kernel))
end
```

```@example hmf
dev = CPU()
nsteps = 1000
dt = 0.1

mass = 1.0
T = 0.1
mesh1 = OneDGrid(dev, 64, -π, π)
mesh2 = OneDGrid(dev, 128, -10, 10)

n1, delta1 = mesh1.len, mesh1.step
n2, delta2 = mesh2.len, mesh2.step
x, v = mesh1.points, mesh2.points'
ϵ = 0.1

b = 1 / T
m = mag(b, mass)

w   = sqrt(m)
f   = zeros(Complex{Float64}, (n1,n2))
fᵗ  = zeros(Complex{Float64}, (n2,n1))

@. f  = exp(-b * ((v^2 / 2) - m * cos(x)))

a   = mass / Norm(real(f), delta1, delta2)
@. f =  a * exp(-b * ((v^2 / 2) - m * cos(x))) * (1 + ϵ * cos(x))

transpose!(fᵗ, f)
@show size(f), size(fᵗ)
contour(mesh1.points, mesh2.points, real(fᵗ))
```

Compute the electrical energy in ``hmf_poisson!`` and display it.

We also analytically compute E(0, x) to assess the correctness of the numerical computation. We have
```math
E(t,x) = -\partial_x\Phi(t,x),
```
where
```math
\Phi(t,x) = \frac{1}{\pi}\int_{[-\pi, \pi]\times \mathbb{R}} \cos(x-y)f(t,y,v)dydv.
```

By trigonometrical rules, we obtain
```math
\Phi(t,x) = \frac{\cos(x)}{\pi}\int_{[-\pi, \pi]\times \mathbb{R}} \cos(y)f(t,y,v)dydv,
```
where we used the oddness of the mapping $y\mapsto sin(y)f(t,y,v)$ to forget about the second integral that should
appear. Using the formula for $f$
```math
f(0, y, v) = a \exp(-b(v^2/2 - m\cos(y))) (1 + \epsilon\cos(y)),
```
we obtain after some manipulations 
```math
\Phi(0,x) = 2a \sqrt{\frac{2\pi}{b}} \cos(x) \left( I_1(bm) + \frac{\epsilon}{2} (I_0(bm)+I_2(bm)) \right).
```
Here $t\mapsto I_\nu(t)$ denotes the modified Bessel function of the first kind, of order $\nu \in \mathbb{N}$. The multiplicative 
coefficient 2 comes from the definition of $I_\nu$:
```math
I_\nu(t) = \frac{1}{\pi}\int_{[0, \pi]} e^{t\cos(y)}\cos(\nu y)dy = \frac{1}{2\pi} \int_{[-\pi, \pi]} e^{t\cos(y)}\cos(\nu y)dy
```

For the parameters chosen, we obtain that 
```math
\Phi(0,x) = \alpha \cos(x),\, E(0, x) = -\partial_x \Phi(0, x) = \alpha \sin(x)
```
where $\alpha = 0.32962331549891355$.

```@example hmf
ex = zeros(Float64,n1)
hmf_poisson!(fᵗ, mesh1, mesh2, ex)
test = copy(f)
T = Float64[]
plot(x, ex, label="E(t,x)", minorgrid=true)
α = 0.32962331549891355
plot!(x, α*sin.(x), markershape=:circle, linewidth=0, label="αsin(x)")
```

```@example hmf
maximum(ex) - α
```

Hence our numerical computation of ``ex`` through ``hmf_poisson!`` gives the expected result.


## Benchmarks

We can sometimes find in the literature that the Vlasov-HMF model is a good simplified version of the Vlasov-Poisson
equations, since it models the expected behavior quite remarkably. Another way of saying it is that after a short time
only the Fourier modes $k=\pm 1$ remain, the others vanishing rapidly. 

To  illustrate this numerically, we apply the above code -- which solves the Vlasov-HMF equations -- to the examples 
used in the Vlasov-Poisson case. The theoretical analysis for the analytical results can be found in the book 
Numerical Methods for the Vlasov-Maxwell Equations from E. Sonnendrücker.

The results are to be compared with those obtained by solving the Vlasov-Poisson equations.



### Landau damping


Here we consider $x\in[0, 2\pi/k_x]$, with $k_x$ some parameter. The initial condition reads
```math
f_0(x,v) = (1+\epsilon \cos(k_x x)) \frac{e^{-v^2/2}}{\sqrt{2\pi}}



#### $k_x = 0.5$

```@example hmf
import VlasovSolvers: advection!

dev = CPU()
nsteps = 500
dt = 0.1

ϵ = 0.01
kx = 0.5
# For this set of parameters the expected damping coefficient is -0.1533,
# and the "period of the damping" is ≈2π/1.4156.

mesh1 = OneDGrid(dev, 64, 0, 2π / kx)
mesh2 = OneDGrid(dev, 128, -10, 10)

n1, delta1 = mesh1.len, mesh1.step
n2, delta2 = mesh2.len, mesh2.step
x, v = mesh1.points, mesh2.points'

flandau  = zeros(Complex{Float64}, (n1,n2))
flandauᵗ = zeros(Complex{Float64}, (n2,n1))

@. flandau  = exp(-v^2 / 2) * (1 + ϵ * cos(x*kx)) / sqrt(2π)
transpose!(flandauᵗ, flandau)

Elandau = Array{Float64}(undef, nsteps)

ex = zeros(Float64,n1)
hmf_poisson!(flandauᵗ, mesh1, mesh2, ex)
advection!(flandauᵗ, mesh2, ex, 0.5dt)


for n in 1:nsteps
    Elandau[n] = sqrt(sum(ex.^2) * mesh1.step)

    advection!(flandau, mesh1, v, dt)
    transpose!(flandauᵗ, flandau)
    hmf_poisson!(flandauᵗ, mesh1, mesh2, ex)
    advection!(flandauᵗ, mesh2, ex, dt)
    transpose!(flandau, flandauᵗ)
end
```

```@example hmf
t = range(0., stop=nsteps*dt, length=nsteps) |> collect

period = 2π / 1.4156
start = 10.5
nbperiods = [2, 4, 6, 8]
vline([start], ls=:dash, label="Period reference", c="black")
for p=nbperiods
  vline!([start + p.*period], ls=:dash, label="After $p periods", c="grey")
end

plot!(t, log.(Elandau), xlabel = "t", minorgrid=true, label="log(E(t))")
plot!(x-> - 0.1533x - 3.2, label="y = -0.1533x - 5.5")
# Check "periods" of energy : theoretically they are ≈2π/1.4156.
# Each period in the energy is materialized by 2 bumps.

plot!(t, log.(abs.(4ϵ.*0.3677.*exp.(-0.1533.*t).*cos.(1.4156.*t .- 0.536245)) * √(mesh1.stop/2.)), label="theoretical", ls=:dash)
plot!(legend=:topright)
```

#### $k_x = 0.4$ 

```@example hmf
import VlasovSolvers: advection!

dev = CPU()
nsteps = 500
dt = 0.1

ϵ = 0.001
kx = 0.4
# For this set of parameters the expected damping coefficient is -0.0661,
# and the "period of the damping" is ≈2π/1.2850.

mesh1 = OneDGrid(dev, 64, 0, 2π / kx)
mesh2 = OneDGrid(dev, 128, -10, 10)

n1, delta1 = mesh1.len, mesh1.step
n2, delta2 = mesh2.len, mesh2.step
x, v = mesh1.points, mesh2.points'

flandau  = zeros(Complex{Float64}, (n1,n2))
flandauᵗ = zeros(Complex{Float64}, (n2,n1))

@. flandau  = exp(-v^2 / 2) * (1 + ϵ * cos(x*kx)) / sqrt(2π)
transpose!(flandauᵗ, flandau)

Elandau = Array{Float64}(undef, nsteps)

ex = zeros(Float64,n1)
hmf_poisson!(flandauᵗ, mesh1, mesh2, ex)
advection!(flandauᵗ, mesh2, ex, 0.5dt)


for n in 1:nsteps
    Elandau[n] = sqrt(sum(ex.^2) * mesh1.step)

    advection!(flandau, mesh1, v, dt)
    transpose!(flandauᵗ, flandau)
    hmf_poisson!(flandauᵗ, mesh1, mesh2, ex)
    advection!(flandauᵗ, mesh2, ex, dt)
    transpose!(flandau, flandauᵗ)
end
```

```@example hmf
t = range(0., stop=nsteps*dt, length=nsteps) |> collect

period = 2π / 1.2850
start = 8.7
nbperiods = [2, 4, 6, 8]
vline([start], ls=:dash, label="Period reference", c="black")
for p=nbperiods
  vline!([start + p.*period], ls=:dash, label="After $p periods", c="grey")
end


plot!(t, log.(Elandau), xlabel = "t", minorgrid=true, label="log(E(t))")
plot!(x-> - 0.0661x - 5.2, label="y = -0.0661x - 5.2")
# Check "periods" of energy : theoretically they are ≈2π/1.4156.
# Each period in the energy is materialized by 2 bumps.

plot!(t, log.(4ϵ.*0.42466.*exp.(-0.0661.*t).*abs.(cos.(1.285.*t .- 0.33577)) .* √((mesh1.stop - mesh1.start)/2.)), label="theoretical", ls=:dash)
plot!(legend=:topright)
```

### Two-stream Instability


Again, we consider $x\in [0, 2\pi/k_x]$. The initial condition reads 
```math
f_0(x,v) = (1+\epsilon\cos(k_x x)) \frac{e^{-(v+v_0)^2/2} + e^{-(v-v_0)^2/2}}{2\sqrt{2\pi}}
```


#### $k_x=0.2, v_0=1.3$

```@example hmf
import VlasovSolvers: advection!

dev = CPU()
nsteps = 1000
dt = 0.1

ϵ = 0.001
kx = 0.2
v0 = 1.3
# For this set of parameters the expected damping coefficient is -0.00104
# and the "period of the damping" is ≈2π/1.1648.

mesh1 = OneDGrid(dev, 64, 0, 2π / kx)
mesh2 = OneDGrid(dev, 128, -10, 10)

n1, delta1 = mesh1.len, mesh1.step
n2, delta2 = mesh2.len, mesh2.step
x, v = mesh1.points, mesh2.points'

ftsi  = zeros(Complex{Float64}, (n1,n2))
ftsiᵗ  = zeros(Complex{Float64}, (n2,n1))

@.  ftsi =  (exp(- 0.5 * (v-v0)^2) + exp(-0.5 * (v+v0)^2)) * (1 + ϵ * cos(kx*x)) / (2*√(2π))
transpose!(ftsiᵗ, ftsi)

Etsi = Array{Float64}(undef, nsteps)

ex = zeros(Float64,n1)
hmf_poisson!(ftsiᵗ, mesh1, mesh2, ex)
advection!(ftsiᵗ, mesh2, ex, 0.5dt)


for n in 1:nsteps
    Etsi[n] = sqrt(sum(ex.^2) * mesh1.step)

    advection!(ftsi, mesh1, v, dt)
    transpose!(ftsiᵗ, ftsi)
    hmf_poisson!(ftsiᵗ, mesh1, mesh2, ex)
    advection!(ftsiᵗ, mesh2, ex, dt)
    transpose!(ftsi, ftsiᵗ)
end
```

```@example hmf
t = range(0., stop=nsteps*dt, length=nsteps) |> collect

plot(t, log.(Etsi), xlabel = "t", minorgrid=true, label="log(E(t))")
plot!(x-> -0.00104x - 4.15, label="y = -0.00104x - 4.15")
# Check "periods" of energy : theoretically they are ≈2π/1.1648.
# Each period in the energy is materialized by 2 bumps.
period = 2π / 1.1648
start = 9.5
nbperiods = [5, 10, 15]
vline!([start], ls=:dash, label="Period reference", c="black")
for p=nbperiods
  vline!([start + p.*period], ls=:dash, label="After $p periods", c="grey")
end

plot!(legend=:bottomright)
```

#### $k_x = 0.2, v_0 = 3.0$

```@example hmf
import VlasovSolvers: advection!

dev = CPU()
nsteps = 500
dt = 0.1

ϵ = 0.001
kx = 0.2
v0 = 3.0
# For this set of parameters the expected damping coefficient is +0.2845.

mesh1 = OneDGrid(dev, 64, 0, 2π / kx)
mesh2 = OneDGrid(dev, 128, -10, 10)

n1, delta1 = mesh1.len, mesh1.step
n2, delta2 = mesh2.len, mesh2.step
x, v = mesh1.points, mesh2.points'

ftsi  = zeros(Complex{Float64}, (n1,n2))
ftsiᵗ  = zeros(Complex{Float64}, (n2,n1))

@.  ftsi =  (exp(- 0.5 * (v-v0)^2) + exp(-0.5 * (v+v0)^2)) * (1 + ϵ * cos(kx*x)) / (2*√(2π))
transpose!(ftsiᵗ, ftsi)

Etsi = Array{Float64}(undef, nsteps)

ex = zeros(Float64,n1)
hmf_poisson!(ftsiᵗ, mesh1, mesh2, ex)
advection!(ftsiᵗ, mesh2, ex, 0.5dt)


for n in 1:nsteps
    Etsi[n] = sqrt(sum(ex.^2) * mesh1.step)

    advection!(ftsi, mesh1, v, dt)
    transpose!(ftsiᵗ, ftsi)
    hmf_poisson!(ftsiᵗ, mesh1, mesh2, ex)
    advection!(ftsiᵗ, mesh2, ex, dt)
    transpose!(ftsi, ftsiᵗ)
end
```

```@example hmf
t = range(0., stop=nsteps*dt, length=nsteps) |> collect

plot(t, log.(Etsi), xlabel = "t", minorgrid=true, label="log(E(t))")
plot!(x-> 0.2845x - 6.2, label="y = 0.2845x - 6.2")
```

## Conclusion


We illustrated the fact that the Vlasov-HMF equations is a toy model that exhibits a behavior very close to the Vlasov-Poisson system. 
