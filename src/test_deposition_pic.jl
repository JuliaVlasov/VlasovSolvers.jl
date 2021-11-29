#=
Pour vérifier la déposition on sample $x\mapsto \cos(x)$, puis on dépose sur la grille, et on résout l'équation de Poisson $-\Delta \Phi = \rho$. On doit retrouver du cosinus.
=#

using LinearAlgebra, QuadGK, Roots, FFTW
using VlasovSolvers
using Plots
using ProgressMeter
import VlasovSolvers: advection!

include("PIC.jl")

function samples_verif1(nsamples, kx, α::Float64)
    #=
    Sample distribution defined by fx in space and a gaussian of variance beta and mean μ in velocity.
    =#
    x0 = Array{Float64}(undef, nsamples)
    seq = SobolSeq(1)
    @showprogress for i in 1:nsamples
        r = Sobol.next!(seq)
        resx = find_zero(x-> x + α/kx * sin(kx*x) - r[1]*2π/kx, 0.0)
        x0[i] = resx
    end
    v0 = rand(Uniform(-10, 10), nsamples)
    wei = 2π/kx/nsamples .* ones(Float64, nsamples)
    return x0, v0, wei
end

nstep = 1000
dt = 0.05
kx = 0.5
L = 2π / kx
K = 1 # paramètre de troncature du noyau
np = 100000 # nb de particules
nx = 32   # nb de mailles pour le calcul de l energie
ϵ = 0.5 # perturbation initiale
μ = 0.0
β = 1.0


@time (x0, y0, wei) = samples_verif1(np, kx, ϵ)
p = Particles(x0, y0, wei, np);

