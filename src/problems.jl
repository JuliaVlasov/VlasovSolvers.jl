export VlasovProblem

abstract type AbstractProblem end

struct VlasovProblem{Method<:AbstractMethod} <: AbstractProblem

    f :: DistributionFunction
    method :: Method
    dev :: AbstractDevice

    function VlasovProblem( f, method:: BSLSpline, dev)

        new{BSLSpline}( f, method, dev)

    end

    function VlasovProblem( f, method::Fourier, dev)

        new{Fourier}( f, method, dev)
      
    end

end

export solve!

function solve!( problem::VlasovProblem{BSLSpline}, stepper::StrangSplitting, dt, nsteps )

  nrj = Float64[]
  
  for it in 1:nsteps
        
     advection_x!( problem.f, 0.5dt)
     sol = advection_v!( problem.f, dt)
     push!(nrj, sol)
     advection_x!( problem.f, 0.5dt)
        
  end
                  
  nrj

end

function solve!( problem::VlasovProblem{Fourier}, stepper :: StrangSplitting, dt, nsteps)

    # Initialize distribution function
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
    
    nrj = Float64[]
    
    for i in 1:nsteps

        advection_v!(fᵀ, problem.method, e,  0.5dt)
        transpose!(f,fᵀ)
        advection_x!( f, problem.method, e, v, dt)
        push!(nrj, log(sqrt((sum(e.^2)) * dx)))
        transpose!(fᵀ,f)
        advection_v!(fᵀ, problem.method, e,  0.5dt)

    end

    nrj

end


