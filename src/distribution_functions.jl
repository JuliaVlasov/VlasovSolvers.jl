using FFTW
using LinearAlgebra
using Statistics

export DistributionFunction

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct DistributionFunction

    xgrid :: OneDGrid
    vgrid :: OneDGrid
    values :: Array{AbstractFloat, 2}

    function DistributionFunction( xgrid, vgrid )

        nx = xgrid.len
        nv = vgrid.len
        values = zeros(nx, nv)

        new( xgrid, vgrid, values )

    end

end

export landau!

"""
$(SIGNATURES)

Initialize the distribution function for the Landau damping test case

```math
f(x,v,t=0) = \\frac{1}{\\sqrt{2\\pi}} ( 1 + \\cos{kx \cdot x} ) \\exp (-\\frac{v^2}{2})
```

"""
function landau!( f :: DistributionFunction, α, kx)

    nx = f.xgrid.len
    nv = f.vgrid.len
    x  = f.xgrid.points
    v  = f.vgrid.points
    f.values .= (1 .+ α*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.^2))

end

export two_stream_instability!

"""
$(SIGNATURES)
"""
function two_stream_instability!(f; eps = 0.01, xi = 0.90, v0 = 2.4)

    nx = f.xgrid.len
    nv = f.vgrid.len
    xg = f.xgrid.points
    vg = f.vgrid.points

    for (i,x) in enumerate(xg), (j,v) in enumerate(vg)
        f.values[i,j] = ((1 + eps*((cos(4π*x) + cos(6π*x)) / 1.2 
                .+ cos(2π*x))) * (1/sqrt(2π)) * ((2-2xi)/(3-2xi))
                    * (1 + 5 * v^2 / (1-xi)) * exp(-.5 * v^2))
    end
end 

