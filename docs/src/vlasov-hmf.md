---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.1
  kernelspec:
    display_name: Julia 1.6.3
    language: julia
    name: julia-1.6
---

In this example we will show how the VlasovSolvers package can be used to solve the Vlasov equations in the Hamiltonian 
Mean Field framework (usually called Vlasov-HMF).

Just like the Vlasov-Poisson equations are the Vlasov-Maxwell equations where field equations are simplified, the 
Vlasov-HMF equations are the Vlasov-Poisson equations where the field equations are further simplified.

Recall the Vlasov equations in its most simplified case (i.e. no source term):
$$
\partial_t f + v\cdot \partial_x f + E\cdot \partial_v f = 0
$$

The quantity $E$ is the electric field, defined in the Poisson framework by
$$
-\Delta \phi = 1 - \rho,\ E = -\nabla \phi,\ \rho(t,x) = \int f(t,x,v)dv
$$

The Poisson equation $-\Delta \phi = 1 - \rho$ on a periodic space-domain $[0, L]$ is usually solved by means of a Fourier
transform. In the discrete case this is performed by a DFT, involving $N_x$ Fourier modes in total (where $N_x$ is the
number of spatial nodes). In the HMF framework, we apply the same idea but restrict ourselves to the Fourier modes 
corresponding to the frequencies $k=\pm 1$ (the mode corresponding to $k=0$ is zero since $\phi$ has a zero average).

```julia
using LinearAlgebra, QuadGK, Roots, FFTW
using VlasovSolvers
using Plots
```

```julia
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

```julia
function Norm(f::Array{Float64,2}, delta1, delta2)
    delta1 * sum(delta2 * sum(real(f), dims=1))
end
```

As explained above, we compute the electric field `ex` in a similar manner as in Poisson case, but the kernel
corresponding to solving the Laplacian in periodic space domain is restricted to its modes $k=\pm 1$.

```julia
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
    rho = mesh2.step .* vec(sum(fᵗ, dims=1)) # ≈ ∫ f(t,x_i,v)dv, i=1, ..., n1
    kernel = zeros(Float64, n1)
    k = -(mesh1.stop - mesh1.start) / (2π)
    kernel[2]   =  k    # fourier mode  1
    kernel[end] = -k    # fourier mode -1
    ex .= real(ifft(fft(rho) .* 1im .* kernel))
end
```

```julia
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

f  .= exp.(-b .* ((v.^2 / 2) .- m .* cos.(x)))
transpose!(fᵗ, f)
@show size(f), size(fᵗ)
contour(mesh1.points, mesh2.points, real(fᵗ))
```

```julia
a   = mass / Norm(real(f), delta1, delta2)

@.  f .=  a * exp(-b * (((v^2) / 2) - m * cos(x))) * (1 + ϵ * cos(x))

ex = zeros(Float64,n1)
hmf_poisson!(fᵗ, mesh1, mesh2, ex )
test = copy(f)
T = Float64[]
plot(x, ex)
```

```julia
import VlasovSolvers: advection!

advection!(fᵗ, mesh2, ex, 0.5dt)

for n in 1:nsteps

    gamma1 = Norm(real(f) .* cos.(x), delta1, delta2)
    push!(T,gamma1)

    advection!(f, mesh1, v, dt)
    transpose!(fᵗ, f)
    hmf_poisson!(fᵗ, mesh1, mesh2, ex)
    advection!(fᵗ, mesh2, ex, dt)
    transpose!(f, fᵗ)

end

# Substracting from gamma its long time average

Gamma1 = Norm(real(f) .* cos.(x), delta1, delta2)

T .= T .- Gamma1

t = range(0., stop=nsteps*dt, length=nsteps) |> collect
T .= abs.(T)

plot(t, log.(T), xlabel = "t", ylabel = "|C[f](t)-C[f][T]|")
```

# Benchmark

We now compare our numerical results with theoretical ones from [Sonnendrücker, Numerical For the Vlasov-Maxwell Equations, pp.54-59]


## Landau damping:


### $k=0.5$

```julia
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

```julia
t = range(0., stop=nsteps*dt, length=nsteps) |> collect

period = 2π / 1.4156
start = 10.5
minElandau, maxElandau = extrema(log.(Elandau))
vline(start .+ [0, 5, 10].*period, ls=:dash)

plot!(t, log.(Elandau), xlabel = "t", minorgrid=true, label="log(energy)")
plot!(x-> - 0.1533x - 3.2, label="y = -0.1533x - 3.2")
# Check "periods" of energy : theoretically they are ≈2π/1.4156. 
# Each period in the energy is materialized by 2 bumps.

plot!(t, log.(4ϵ.*0.3677.*exp.(-0.1533.*t).*abs.(cos.(1.4156.*t .- 0.536245)) * √((mesh1.stop - mesh1.start)/2.)), 
    label="log(4ϵ*0.3677*exp(-0.1533t)*|cos(1.4156t - 0.536245)| * √(L/2)")
plot!(legend=:topright)
```

### $k=0.4$

```julia
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

```julia
t = range(0., stop=nsteps*dt, length=nsteps) |> collect

period = 2π / 1.2850
start = 8.7
minElandau, maxElandau = extrema(log.(Elandau))
vline(start .+ [0, 5].*period, ls=:dash)

plot!(t, log.(Elandau), xlabel = "t", minorgrid=true, label="log(energy)")
plot!(x-> - 0.0661x - 5.2, label="y = -0.0661x - 5.2")
# Check "periods" of energy : theoretically they are ≈2π/1.4156. 
# Each period in the energy is materialized by 2 bumps.

plot!(t, log.(4ϵ.*0.42466.*exp.(-0.0661.*t).*abs.(cos.(1.285.*t .- 0.33577)) .* √((mesh1.stop - mesh1.start)/2.)),
  label="log(4ϵ*0.42466*exp(-0.0661t) |(cos(1.285t - 0.33577)| * √(L/2))")
plot!(legend=:topright)
```

## Two-stream instability:


### $k=0.2,\ v_0 = 1.3$

```julia
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

```julia
t = range(0., stop=nsteps*dt, length=nsteps) |> collect

plot(t, log.(Etsi), xlabel = "t", minorgrid=true, label="log(energy)")
plot!(x-> -0.00104x - 4.15, label="y = -0.00104x - 4.15")
# Check "periods" of energy : theoretically they are ≈2π/1.1648. 
# Each period in the energy is materialized by 2 bumps.
period = 2π / 1.1648
start = 9.5
minEtsi, maxEtsi = extrema(log.(Etsi))
plot!([start, start], [minEtsi, maxEtsi], label="period reference", ls=:dash)
plot!([start+5period, start+5period], [minEtsi, maxEtsi], label="after 5 theoretical periods", ls=:dash) # Should cover exactly 10 bumps
plot!([start+10period, start+10period], [minEtsi, maxEtsi], label="after 10 theoretical periods", ls=:dash) # Should cover exactly 20 bumps
plot!(legend=:bottomright)
```

### $k=0.2,\ v_0 = 3.0$

```julia
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

```julia
t = range(0., stop=nsteps*dt, length=nsteps) |> collect

plot(t, log.(Etsi), xlabel = "t", minorgrid=true, label="log(energy)")
plot!(x-> 0.2845x - 6.2, label="y = 0.2845x - 6.2")
```
