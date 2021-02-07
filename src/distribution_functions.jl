using FFTW
using LinearAlgebra
using Statistics

export DistributionFunction

"""
    DistributionFunction( grid1, grid2 )

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
    landau!( f, α, kx )

Initialize the distribution function for the Landau damping test case

```math
f(x,v,t=0) = \\frac{1}{\\sqrt{2\\pi}} ( 1 + \\cos (kx x) ) \\exp (-\\frac{v^2}{2})
```

"""
function landau!( f :: DistributionFunction, α, kx)

    nx = f.xgrid.len
    nv = f.vgrid.len
    x  = f.xgrid.points
    v  = f.vgrid.points
    f.values .= (1 .+ α*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.^2))

end



"""
    bspline(p, j, x)

Return the value at x in [0,1[ of the B-spline with integer nodes of degree p with support starting at j.
Implemented recursively using the [De Boor's Algorithm](https://en.wikipedia.org/wiki/De_Boor%27s_algorithm)

```math
B_{i,0}(x) := \\left\\{
\\begin{matrix}
1 & \\mathrm{if}  \\quad t_i ≤ x < t_{i+1} \\\\
0 & \\mathrm{otherwise} 
\\end{matrix}
\\right.
```

```math
B_{i,p}(x) := \\frac{x - t_i}{t_{i+p} - t_i} B_{i,p-1}(x) 
+ \\frac{t_{i+p+1} - x}{t_{i+p+1} - t_{i+1}} B_{i+1,p-1}(x).
```
"""
function bspline(p::Int, j::Int, x::Float64)
   if p == 0
       if j == 0
           return 1.0
       else
           return 0.0
       end
   else
       w = (x - j) / p
       w1 = (x - j - 1) / p
   end
   return (w * bspline(p - 1, j, x) + (1 - w1) * bspline(p - 1, j + 1, x))
end


"""
    advection!(f, grid, v, dt; p = 5)

function to advect the distribution function `f` with velocity `v` along 
first `f` dimension during a time step `dt`. Interpolation method uses bspline periodic of order 5 by default. Real type version.
"""
function advection!(f    :: Array{AbstractFloat, 2},
                    grid :: OneDGrid, 
                    v, 
                    dt;
                    p = 5)
    
   nx = grid.len
   nv = length(v)
   dx = grid.step
   modes = [2π * i / nx for i in 0:nx-1]
    
   # compute eigenvalues of degree p b-spline matrix
   eig_bspl  = zeros(Float64, nx)
   eig_bspl .= bspline(p, -div(p+1,2), 0.0)
   for i in 1:div(p+1,2)-1
      eig_bspl .+= bspline(p, i - (p+1)÷2, 0.0) * 2 .* cos.(i * modes)
   end
   eigalpha = zeros(Complex{Float64}, nx)
    
   ft = fft(f,1)
    
   for j in 1:nv
      alpha = dt * v[j] / dx
      # compute eigenvalues of cubic splines evaluated 
      # at displaced points
      ishift = floor(-alpha)
      beta   = -ishift - alpha
      fill!(eigalpha,0.0im)
      for i in -div(p-1,2):div(p+1,2)
         eigalpha .+= (bspline(p, i-div(p+1,2), beta) 
                        .* exp.((ishift+i) * 1im .* modes))
      end
          
      # compute interpolating spline using fft and properties 
      # of circulant matrices
      
      ft[:,j] .*= eigalpha ./ eig_bspl
        
   end
        
   f .= real(ifft(ft,1))
    
end            

"""
    advection!(f, grid, v, dt; p = 5)

function to advect the distribution function `f` with velocity `v` along 
first `f` dimension during a time step `dt`. Interpolation method uses bspline periodic of order 5 by default. Complex type version.
"""
function advection!(f    :: Array{ComplexF64, 2},
                    grid :: OneDGrid, 
                    v, 
                    dt;
                    p = 5)
    
   nx = grid.len
   nv = length(v)
   dx = grid.step
   modes = [2π * i / nx for i in 0:nx-1]
    
   # compute eigenvalues of degree p b-spline matrix
   eig_bspl  = zeros(Float64, nx)
   eig_bspl .= bspline(p, -div(p+1,2), 0.0)
   for i in 1:div(p+1,2)-1
      eig_bspl .+= bspline(p, i - (p+1)÷2, 0.0) * 2 .* cos.(i * modes)
   end
   eigalpha = zeros(Complex{Float64}, nx)
    
   fft!(f,1)
    
   for j in 1:nv
      alpha = dt * v[j] / dx
      # compute eigenvalues of cubic splines evaluated 
      # at displaced points
      ishift = floor(-alpha)
      beta   = -ishift - alpha
      fill!(eigalpha,0.0im)
      for i in -div(p-1,2):div(p+1,2)
         eigalpha .+= (bspline(p, i-div(p+1,2), beta) 
                        .* exp.((ishift+i) * 1im .* modes))
      end
          
      # compute interpolating spline using fft and properties 
      # of circulant matrices
      
      f[:,j] .*= eigalpha ./ eig_bspl
        
   end
        
   ifft!(f,1)
    
end            


"""
    compute_rho(f)

Compute charge density
ρ(x,t) = ∫ f(x,v,t) dv
"""
function compute_rho( f::DistributionFunction)
    
   dv = f.vgrid.step
   rho = dv * sum(real(f.values), dims=2)
   vec(rho .- mean(rho)) # vec squeezes the 2d array returned by sum function

end

# + {"slideshow": {"slide_type": "slide"}}
"""
    compute_e(f)

compute Ex using that -ik*Ex = rho 
"""
function compute_e( f::DistributionFunction )

   rho = compute_rho( f )
   nx = f.xgrid.len
   xmin = f.xgrid.start
   xmax = f.xgrid.stop
   k =  2π / (xmax - xmin)
   modes = zeros(Float64, nx)
   modes .= k * vcat(0:div(nx,2)-1,-div(nx,2):-1)
   modes[1] = 1.0
   rhok = fft(rho)./modes
   rhok .*= -1im
   ifft!(rhok)
   real(rhok)

end


function advection_x!( f :: DistributionFunction, dt) 

    advection!(f.values, f.xgrid, collect(f.vgrid.points), dt)

end

function advection_v!( f :: DistributionFunction, dt)

    dx = f.xgrid.step
    nx = f.xgrid.len
    e  = compute_e(f)
    nrj = 0.5*log(sum(e.*e)*dx)
    fᵗ = collect(transpose(f.values))
    advection!(fᵗ, f.vgrid, e, dt)
    f.values .= transpose(fᵗ)
    nrj

end

