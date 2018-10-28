#!/usr/bin/env julia

using Pkg

Pkg.add("https://github.com/JuliaVlasov/VlasovBase.jl")
Pkg.add("https://github.com/JuliaVlasov/SemiLagrangian.jl")
Pkg.add("https://github.com/JuliaVlasov/Fourier.jl")
Pkg.add("https://github.com/JuliaVlasov/SplittingOperators.jl")
Pkg.add("https://github.com/JuliaVlasov/VlasovExamples.jl")

using VlasovSolvers, Test

@test true
