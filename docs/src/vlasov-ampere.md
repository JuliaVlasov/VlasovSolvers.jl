```@meta
EditURL = "<unknown>/docs/examples/vlasov-ampere.jl"
```

# Vlasov-Ampere

[`notebook`](<unknown>notebooks/vlasov-ampere.ipynb)

Compute Landau damping by solving Vlasov-Ampere system.

```math
 \frac{∂f}{∂t} + υ \frac{∂f}{∂x} - E(t,x) \frac{∂f}{∂υ} = 0
```

```math
\frac{∂E}{∂t} = - J = ∫ fυ dυ
```

## Algorithm

 - For each ``j`` compute discrete Fourier transform in ``x`` of
   ``(x_i,υ_j)`` yielding ``f_k^n(υ_j)``,
 - For `` k ≂̸ 0 ``
     * Compute ``f^{n+1}_k(υ_j) = e^{−2iπ k υ Δt/L} f_n^k(υ_j), ``
     * Compute ``ρ_k^{n+1} = Δ υ ∑_j􏰄 f^{n+1}_k(υ_j), ``
     * Compute ``E^{n+1}_k = ρ^{n+1}_k L/(2iπkϵ_0), ``
 - For ``k = 0`` do nothing:
```math
f_{n+1}(υ_j) = f^n_k(υ_j), E^{n+1}_k = E^n_k.
```
 - Perform in2erse discrete Fourier transform of ``E^{n+1}_k`` and for each
   ``j`` of ``f^{n+1}_k (υ_j)``.

```@example vlasov-ampere
import SemiLagrangian: advection!
import Fourier: Ampere
import SemiLagrangian: advection!,
import SplittingOperators: @Strang
import VlasovBase: compute_rho, compute_e, UniformMesh
using Plots, LinearAlgebra
pyplot()
```

```@example vlasov-ampere
function push_t!( f, fᵀ, mesh1, mesh2, e,  dt)

    advection!( f, fᵀ, mesh1, mesh2, e,  dt, Ampere(), 1 )

end
```

```@example vlasov-ampere
function push_v!(f, fᵀ, mesh1, mesh2, e, dt)

    advection!( f, fᵀ, mesh1, mesh2, e, dt, Ampere(), 2)

end
```

```@example vlasov-ampere
function vm1d( n1, n2, x1min, x1max, x2min, x2max , tf, nt)

    mesh1 = UniformMesh(x1min, x1max, n1, endpoint=false)
    mesh2 = UniformMesh(x2min, x2max, n2, endpoint=false)

    x = mesh1.points
    v = mesh2.points
    ϵ, kx = 0.001, 0.5

    f = zeros(Complex{Float64},(n1,n2))
    fᵀ= zeros(Complex{Float64},(n2,n1))

    f .= (1.0.+ϵ*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.*v))
    transpose!(fᵀ,f)

    e = zeros(Complex{Float64},n1)

    ρ  = compute_rho(mesh2, f)
    e .= compute_e(mesh1, ρ)

    nrj = Float64[]

    dt = tf / nt

    for i in 1:nt

	push!(nrj, 0.5*log(sum(real(e).^2)*mesh1.step))

        @Strang(  push_v!(f, fᵀ, mesh1, mesh2, e, dt),
                  push_t!(f, fᵀ, mesh1, mesh2, e, dt)
               )

    end
    nrj
end
```

```@example vlasov-ampere
n1, n2 = 32, 64
x1min, x1max =  0., 4π
x2min, x2max = -6., 6.
tf = 50
nt = 500

t = range(0,stop=tf,length=nt)
nrj = vm1d(n1, n2, x1min, x1max, x2min, x2max, tf, nt)
plot(t, nrj)
plot!(t, -0.1533*t.-5.50)
savefig("va-plot.png"); nothing
```

![](va-plot.png)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

