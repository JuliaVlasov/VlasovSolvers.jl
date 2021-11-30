#=
We check the deposition step (i.e. computing \rho) works as expected.
=#
using Test

@testset "PIC deposition : uniform    " begin
    using LinearAlgebra
    using VlasovSolvers
    # using Plots

    include("PIC.jl")

    function samples_verif1(nsamples, kx::Float64)
        x0 = Array{Float64}(undef, nsamples)
        seq = SobolSeq(1)
        for i = 1:nsamples
            r = Sobol.next!(seq)
            x0[i] = 2π/kx * r[1]
        end
        v0 = rand(Uniform(-10, 10), nsamples)
        wei = 2π/kx/nsamples .* ones(Float64, nsamples)
        return x0, v0, wei
    end

    kx = 0.5
    L = 2π / kx
    np = 100000 # nb de particules
    nx = 32   # nb de mailles pour le calcul de l energie
    mesh1 = OneDGrid(CPU(), nx, 0, L)

    (x0, v0, w0) = samples_verif1(np, kx)
    p = Particles(x0, v0, w0, np)

    ρ, ρ_tot = compute_rho(p, mesh1)

    # arrx = mesh1.start:mesh1.step:(mesh1.stop - mesh1.step)
    # plot(arrx, ρ)

    @test maximum(abs.(ρ .- 1.)) < 10^-3
    @test ρ_tot ≈ 1.0
end

@testset "PIC deposition : 1 + 0.99cos" begin
    using LinearAlgebra, Roots
    using VlasovSolvers
    # using Plots

    include("PIC.jl")

    function samples_verif2(nsamples)
        x0 = Array{Float64}(undef, nsamples)
        seq = SobolSeq(1)
        for i = 1:nsamples
            r = Sobol.next!(seq)
            resx = find_zero(x-> x + 0.99 / kx * sin(kx*x) - r[1]*2π/kx, 0.0)
            x0[i] = resx
        end
        v0 = rand(Uniform(-10, 10), nsamples)
        wei = 2π/kx/nsamples .* ones(Float64, nsamples)
        return x0, v0, wei
    end

    kx = 0.5
    L = 2π / kx
    np = 100000 # nb de particules
    nx = 128   # nb de mailles pour le calcul de l energie
    mesh1 = OneDGrid(CPU(), nx, 0, L)

    (x0, v0, w0) = samples_verif2(np)
    p = Particles(x0, v0, w0, np)

    ρ, ρ_tot = compute_rho(p, mesh1)

    arrx = mesh1.start:mesh1.step:(mesh1.stop - mesh1.step)
    
    @test maximum(abs.(ρ - (1 .+ 0.99cos.(kx.*arrx)))) < 5.10^-3
    @test ρ_tot ≈ 1.0
end