"""
$(SIGNATURES)

Advect the distribution function `f` with velocity `v` along first `f` dimension
with a time step `dt`. Interpolation method uses bspline periodic of order 5 by default. Real type version.
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
      fill!(eigalpha, 0.0im)
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
$(SIGNATURES)

Advect the distribution function `f` with velocity `v` along first `f` dimension
with a time step `dt`. Interpolation method uses bspline periodic of order 5 by default. Complex type version.
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
$(SIGNATURES)

Compute charge density
ρ(x,t) = ∫ f(x,v,t) dv
"""
function compute_rho( f::DistributionFunction)
    
   dv = f.vgrid.step
   rho = dv * sum(real(f.values), dims=2)
   vec(rho .- mean(rho)) # vec squeezes the 2d array returned by sum function

end

"""
$(SIGNATURES)

compute Ex using that -ik*Ex = rho 
"""
function compute_e( f::DistributionFunction )

   rho = compute_rho( f )
   nx = f.xgrid.len
   xmin = f.xgrid.start
   xmax = f.xgrid.stop
   k =  2π / (xmax - xmin)
   modes = k * vcat(0:div(nx,2)-1, -div(nx,2):-1)
   modes[1] = 1.0
   rhok = fft(rho) ./ modes .* (-1im)
   real(ifft!(rhok))

end


function advection_x!( f :: DistributionFunction, dt, method::BSLSpline)

    advection!(f.values, f.xgrid, f.vgrid.points, dt; method.p)

end

function advection_v!( f :: DistributionFunction, dt, method::BSLSpline)

    dx = f.xgrid.step
    nx = f.xgrid.len
    e  = compute_e(f)
    nrj = 0.5*log(sum(e.*e)*dx)
    fᵗ = collect(transpose(f.values))
    advection!(fᵗ, f.vgrid, e, dt; method.p)
    f.values .= transpose(fᵗ)
    nrj

end

function advection_x!( f :: DistributionFunction, dt, method::BSLLagrange)

    nx = f.xgrid.len
    fp = zeros(nx)
    fi = zeros(nx)

    for j in eachindex(f.vgrid.points)

        v = f.vgrid.points[j]
        fi .= f.values[:,j]
        advection!(fp, fi, f.xgrid, method.interp, v, dt)
        f.values[:,j] .= fp

    end

end

function advection_v!( f :: DistributionFunction, dt, method::BSLLagrange)

    dx = f.xgrid.step
    nx = f.xgrid.len
    e  = compute_e(f)
    nrj = 0.5*log(sum(e.*e)*dx)
    
    nv = f.vgrid.len
    fi = zeros(nv)
    fp = zeros(nv)
    for i in eachindex(e)
        fi .= f.values[i,:]
        advection!(fp, fi, f.vgrid, method.interp, e[i], dt)
        f.values[i,:] .= fp
    end

    nrj

end

