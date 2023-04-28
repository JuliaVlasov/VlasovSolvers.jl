export VlasovProblem

abstract type AbstractProblem end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct VlasovProblem{Method<:AbstractMethod} <: AbstractProblem

    f :: DistributionFunction
    method :: Method
    dev :: AbstractDevice

    VlasovProblem( f, method::BSLSpline, dev) = new{BSLSpline}( f, method, dev)

    VlasovProblem( f, method::Fourier, dev) = new{Fourier}( f, method, dev)

    VlasovProblem( f, method::BSLLagrange, dev) = new{BSLLagrange}( f, method, dev)
      
end


"""
$(SIGNATURES)
"""
function solve( problem::AbstractProblem, stepper::StrangSplitting, dt, nsteps)
  
  nrj = Vector{Float64}(undef, nsteps)

  sol = VlasovSolution1D1V() 

  time = 0.0
  for it in 1:nsteps

     advection_x!( problem.f, 0.5dt, problem.method)
     energy = advection_v!( problem.f, dt, problem.method)
     time += dt
     push!(sol.times, time)
     push!(sol.energy, energy)
     advection_x!( problem.f, 0.5dt, problem.method)

  end
                  
  sol

end

"""
$(SIGNATURES)
"""
function solve( problem::VlasovProblem{Fourier}, stepper::StrangSplitting, dt, nsteps)

    x = problem.f.xgrid.points |> collect
    v = problem.f.vgrid.points |> collect
    nx = problem.f.xgrid.len
    nv = problem.f.vgrid.len
    
    f = zeros(Complex{Float64},(nx,nv))
    fᵀ= zeros(Complex{Float64},(nv,nx))
    
    f .= problem.f.values
    
    transpose!(fᵀ,f)
    
    dx  = problem.f.xgrid.step
    dv  = problem.f.vgrid.step
    ρ   = dv .* vec(sum(real(fᵀ), dims=1))
    ρ  .-= mean(ρ)
    e   = zeros(ComplexF64, nx)
    nx  = problem.f.xgrid.len
    Lx  = problem.f.xgrid.stop - problem.f.xgrid.start
    modes  = zeros(Float64, nx)
    modes .= 2π / Lx .* vcat(0:nx÷2-1,-nx÷2:-1)
    modes[1] = 1.0
    ρ̂ = fft(ρ)./modes
    e .= vec(real(ifft(-1im .* ρ̂)))
    
    sol = VlasovSolution1D1V() 
    time = 0.0
    push!(sol.times, time)
    energy = log(sqrt((sum(e.^2)) * dx))
    push!(sol.energy, energy)
    
    for it in 1:nsteps

        advection_v!(fᵀ, problem.method, e, 0.5dt)
        transpose!(f,fᵀ)
        advection_x!( f, problem.method, e, v, dt)
        energy = log(sqrt((sum(e.^2)) * dx))
        time += dt
        push!(sol.times, time)
        push!(sol.energy, energy)
        transpose!(fᵀ,f)
        advection_v!(fᵀ, problem.method, e, 0.5dt)

    end

    sol

end


