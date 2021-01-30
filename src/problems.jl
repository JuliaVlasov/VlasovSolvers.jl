export VlasovProblem

struct VlasovProblem

    f :: DistributionFunction
    method :: AbstractMethod
    dev :: AbstractDevice

end

export solve!

function solve!( problem::VlasovProblem, stepper::StrangSplitting, dt, nsteps )

  nrj = Float64[]
  
  for it in 1:nsteps
        
     advection_x!( problem.f, 0.5dt)
     sol = advection_v!( problem.f, dt)
     push!(nrj, sol)
     advection_x!( problem.f, 0.5dt)
        
  end
                  
  nrj

end

