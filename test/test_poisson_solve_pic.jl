#=
Check that the Poisson solver in the PIC framework works as expected.
It is used in ``PIC.compute_phi``
=#

using Test

@testset "Poisson solve cos" begin
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
    mesh1 = OneDGrid(CPU(), nx, 0, L)

    ρ = cos.(kx .* mesh1.points)
    ρ_tot = 0
    ϕ = compute_phi(ρ, ρ_tot, mesh1)

    # plot(mesh1.points, ρ, label="ρ")
    # plot!(mesh1.points, ϕ .* kx^2, label="ϕ")
    @test maximum(abs.(ϕ  .- ρ ./ kx^2)) < 10^-3
end

@testset "Poisson solve sin" begin
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

    nx = 256
    kx = 0.4
    L = 2π / kx
    mesh1 = OneDGrid(CPU(), nx, 0, L)

    ρ = sin.(5*kx .* mesh1.points)
    ρ_tot = 0
    ϕ = compute_phi(ρ, ρ_tot, mesh1)

    # plot(mesh1.points, ρ, label="ρ")
    # plot!(mesh1.points, ϕ .* 25kx^2, label="ϕ")
    @test maximum(abs.(ϕ  .- ρ ./ (25*kx^2))) < 10^-3
end