```@meta
EditURL = "<unknown>/docs/examples/vlasov-poisson.jl"
```

# Vlasov-Poisson

[notebook](<unknown>notebooks/vlasov-poisson.ipynb)

We consider the dimensionless Vlasov-Poisson equation for one species
with a neutralizing background.

```math
 \frac{∂f}{∂t}+ v⋅∇_x f + E(t,x) ⋅ ∇_v f = 0, \\
 - Δϕ = 1 - ρ, E = - ∇ ϕ \\
 ρ(t,x)  =  ∫ f(t,x,v) dv
```

 - [Vlasov Equation - Wikipedia](https://en.wikipedia.org/wiki/Vlasov_equation)
 - [Landau damping - Wikipedia](https://en.wikipedia.org/wiki/Landau_damping)

```@example vlasov-poisson
using Plots, LinearAlgebra
pyplot()
using VlasovBase
using SplittingOperators
using SemiLagrangian
```

```@example vlasov-poisson
function push_t!(f, mesh1, v, n2, dt)
    advection!(f, mesh1, v, n2, dt, BSpline(5))
end
```

```@example vlasov-poisson
function push_v!(f, fᵗ, mesh1, mesh2, nrj, dt)
    rho = compute_rho(mesh2, f)
    e   = compute_e(mesh1, rho)
    push!(nrj, 0.5*log(sum(e.*e)*mesh1.step))
    transpose!(fᵗ, f)
    advection!(fᵗ, mesh2, e, mesh1.length, dt, BSpline(5))
    transpose!(f, fᵗ)
end
```

```@example vlasov-poisson
function landau(tf, nt)

  n1, n2 = 32, 64
  x1min, x1max = 0.0, 4π
  x2min, x2max = -6., 6.
  mesh1 = UniformMesh(x1min, x1max, n1; endpoint=false)
  mesh2 = UniformMesh(x2min, x2max, n2; endpoint=false)
  x = mesh1.points
  v = mesh2.points

  ϵ, kx = 0.001, 0.5
  f = zeros(Complex{Float64},(n1,n2))
  f .= (1.0.+ϵ*cos.(kx*x))/sqrt(2π) * transpose(exp.(-0.5*v.^2))
  fᵗ = zeros(Complex{Float64},(n2,n1))

  dt = tf / nt

  nrj = Float64[]

  for it in 1:nt
      @Strang( push_t!(f, mesh1, v, n2, dt),
               push_v!(f, fᵗ, mesh1, mesh2, nrj, dt))
  end

  nrj

end
```

```@example vlasov-poisson
nt = 500
tf = 50.0
t  = range(0.0, stop=tf, length=nt)
@time nrj = landau(tf, nt)
plot( t, nrj)
plot!(t, -0.1533*t.-5.50)
savefig("landau-plot.png"); nothing # hide
```

![](landau-plot.png)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

