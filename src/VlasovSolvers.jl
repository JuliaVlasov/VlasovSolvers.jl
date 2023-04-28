module VlasovSolvers

  using DocStringExtensions
  using FFTW
  using LinearAlgebra
  using Statistics

  export solve

  include("devices.jl")
  include("grids.jl")
  include("distribution_functions.jl")
  include("bspline.jl")
  include("lagrange.jl")
  include("advection.jl")
  include("steppers.jl")
  include("fourier.jl")
  include("solution.jl")
  include("problems.jl")

end
