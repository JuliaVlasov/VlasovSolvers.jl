# Vlasov-HMF


```@example 4
using LinearAlgebra, QuadGK, Roots, FFTW
using VlasovSolvers
using Plots
```

```@example 4
" Compute M₀ by solving F(m) = 0 "
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

```@example 4
function Norm(f::Array{Float64,2}, delta1, delta2)
   return delta1 * sum(delta2 * sum(real(f), dims=1))
end
```

```@example 4
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
    rho = mesh2.step .* vec(sum(fᵗ, dims=1))
    kernel = zeros(Float64, n1)
    k = π / (mesh1.stop - mesh1.start)
    kernel[2] = k
    ex .= real(ifft(1im * fft(rho) .* kernel * 4π ))

end
```

```@example 4
dev = CPU()
nsteps = 10000
dt = 0.1

mass = 1.0
T = 0.1
mesh1 = OneDGrid(dev, 64, -π, π)
mesh2 = OneDGrid(dev, 64, -8, 8)

n1, delta1 = mesh1.len, mesh1.step
n2, delta2 = mesh2.len, mesh2.step
x, v = mesh1.points, transpose(mesh2.points)
ϵ = 0.1

b = 1 / T
m = mag(b, mass)

w   = sqrt(m)
f   = zeros(Complex{Float64}, (n1,n2))
fᵗ  = zeros(Complex{Float64}, (n2,n1))

f  .= exp.(-b .* ((v.^2 / 2) .- m .* cos.(x)))
```

```@example 4
a   = mass / Norm(real(f), delta1, delta2)
@.  f =  a * exp(-b * (((v^2) / 2) - m * cos(x))) * (1 + ϵ * cos(x))

ex = zeros(Float64,n1)
hmf_poisson!(f, mesh1, mesh2, ex )
test = copy(f)
T = Float64[]

for n in 1:nsteps

    gamma1 = Norm(real(f) .* cos.(x), delta1, delta2)
    push!(T,gamma1)

    ### Here compute f...

end

#Substracting from gamma its long time average

Gamma1 = Norm(real(f) .* cos.(x), delta1, delta2)

T .= T .- Gamma1

t = range(0., stop=nsteps*dt, length=nsteps) |> collect
T .= abs.(T)

plot(t, log.(T), xlabel = "t", ylabel = "|C[f](t)-C[f][T]|")
```
