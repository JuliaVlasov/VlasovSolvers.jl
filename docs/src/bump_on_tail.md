```@meta
EditURL = "<unknown>/docs/examples/bump_on_tail.jl"
```

# Bump On Tail

[`notebook`](<unknown>notebooks/bump_on_tail.ipynb)

```@example bump_on_tail
import VlasovBase: UniformMesh
import SemiLagrangian: advection!
import VlasovBase: compute_rho, compute_e
import SemiLagrangian: CubicSpline
import SplittingOperators: @Strang

using Plots
using LaTeXStrings

pyplot()
```

```@example bump_on_tail
function push_t!( f, mesh1, v, dt )

    advection!( f, mesh1, v, 0.5dt, CubicSpline(), 1)

end

function push_v!( f, mesh1, mesh2, nrj, dt )

    rho = compute_rho(mesh2, f)
    e   = compute_e(mesh1, rho)
    advection!( f, mesh2, e, dt, CubicSpline(), 2)
    push!(nrj, 0.5*log(sum(e.*e)*mesh1.step))

end

function vlasov_poisson(mesh1  :: UniformMesh,
                        mesh2  :: UniformMesh,
                        f      :: Array{Float64,2},
                        nstep  :: Int64,
                        dt     :: Float64)

    x = mesh1.points
    v = mesh2.points

    nrj = Float64[]
    for istep in 1:nstep

        @Strang(  push_t!( f, mesh1, v, dt ),
                  push_v!( f, mesh1, mesh2, nrj, dt ))


    end
    nrj

end
```

```@example bump_on_tail
α = 0.03
kx  = 0.3
x1min, x1max = 0.0, 2π / kx
n1, n2 = 32, 64
x2min, x2max = -9., 9.
mesh1 = UniformMesh(x1min, x1max, n1)
mesh2 = UniformMesh(x2min, x2max, n2)
f = zeros(Float64,(mesh1.length,mesh2.length))
for (i,x) in enumerate(mesh1.points), (j,v) in enumerate(mesh2.points)
    f[i,j]  = (1.0+α*cos(kx*x)) / (10*sqrt(2π)) * (9*exp(-0.5*v^2)+2*exp(-2*(v-4.5)^2))
end
```

```@example bump_on_tail
nstep = 500
t = range(0.0, stop=50.0, length=nstep)
dt = t[2]
@elapsed nrj = vlasov_poisson( mesh1, mesh2, f, nstep, dt)
```

```@example bump_on_tail
plot(t, nrj, label=L"\frac{1}{2} \log(∫e²dx)")
savefig("bot-plot.png"); nothing # hide
```

![](bot-plot.png)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

