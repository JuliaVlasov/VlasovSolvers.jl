using Sobol
using Roots
# using ProgressMeter
using Random
using Distributions
using SparseArrays
using LinearAlgebra


"""
Describe meta particules, represented by a Dirac disbtibution in (x, v), with a weight (~a high) of wei
"""
struct Particles
    x ::Vector{Float64}     # list of the positions
    v ::Vector{Float64}     # list of the velocities
    wei ::Vector{Float64}   # list of the weights of the particules
    nbpart ::Int64          # nmber of particules
end

"""
Holds pre-allocated arrays
"""
struct ParticleMover
    phi ::Vector{Float64}
    ∂phi :: Vector{Float64}
    meshx
    kx
    K
    C ::Vector{Float64}
    S ::Vector{Float64}
    tmpcosk ::Vector{Float64}
    tmpsink ::Vector{Float64}
    poisson_matrix
    ρ ::Vector{Float64}
    phi_grid ::Vector{Float64}
    
    function ParticleMover(particles::Particles, meshx, kx, K)
        ∂phi = similar(particles.x)
        phi = similar(particles.x)
        tmpcosk = similar(particles.x)
        tmpsink = similar(particles.x)
        nx = meshx.len

        matrix_poisson = spdiagm(  -1 => .- ones(Float64,nx-1),
                                    0 => 2 .* ones(Float64,nx),
                                    1 => .- ones(Float64,nx-1))
        matrix_poisson[1, nx] = -1
        matrix_poisson[nx, 1] = -1

        matrix_poisson ./= meshx.step^2


        new(phi, ∂phi, meshx, kx, K, Vector{Float64}(undef, K), Vector{Float64}(undef, K), tmpcosk, tmpsink
        , matrix_poisson, Vector{Float64}(undef, nx), Vector{Float64}(undef, nx))
    end
end


function samples(nsamples, kx, α::Float64, μ::Float64, β::Float64)
    #=
    Sample distribution defined by 1+α\cos(kx * x) in space and a gaussian of variance beta and mean μ in velocity.
    =#
    x0 = Array{Float64}(undef, nsamples)
    seq = SobolSeq(1)
    for i = 1:nsamples
        r = Sobol.next!(seq)
        resx = find_zero(x-> x + α/kx * sin(kx*x) - r[1]*2π/kx, 0.0)
        x0[i] = resx
    end
    v0 = rand(Normal(μ, 1/√β), nsamples)
    wei = (2π/kx)/nsamples .* ones(Float64, nsamples)
    return x0, v0, wei
end

function update_positions!(p, mesh, dt)
    @inbounds @simd for i = eachindex(p.x)
        p.x[i] += p.v[i] * dt 
        # Periodic boundary conditions
        if p.x[i] >= mesh.stop
            p.x[i] -= mesh.stop
        elseif p.x[i] < mesh.start
            p.x[i] += mesh.stop
        end
    end
end

"""
update particle velocities vp (phi_v)
"""
function update_velocities!(p, pmover, dt)
    compute_S_C!(p, pmover)
    pmover.phi .= 0
    pmover.∂phi .= 0

    @inbounds @simd for k = 1:pmover.K
        pmover.tmpcosk .= cos.(k * pmover.kx .* p.x)
        pmover.tmpsink .= sin.(k * pmover.kx .* p.x)
        denom_phi = pmover.meshx.stop / (2*π^2*k^2) 
        denom_derphi = 1 / (π*k)
        pmover.phi      .+= denom_phi    .* (  pmover.tmpcosk .* pmover.C[k] .+ pmover.tmpsink .* pmover.S[k])
        pmover.∂phi  .+= denom_derphi .* (.-pmover.tmpsink .* pmover.C[k] .+ pmover.tmpcosk .* pmover.S[k])
    end
    p.v .-= dt .* pmover.∂phi
end

""" S[k] = \\sum_l=1^n {\\beta_l * sin(k kx x_l)} et C[k] = \\sum_l=1^n {\\beta_l * cos(k kx x_l)}
utile pour le calcul de la vitesse et du potentiel electrique phi"""
function compute_S_C!(p, pmover)
    pmover.S .= 0
    pmover.C .= 0
    @inbounds @simd for k = 1:pmover.K
        pmover.S[k] = sum(p.wei .* sin.(k .* 2π / pmover.meshx.stop .* p.x))
        pmover.C[k] = sum(p.wei .* cos.(k .* 2π / pmover.meshx.stop .* p.x))
    end
end


function PIC_step!(p::Particles, pmover::ParticleMover, dt)
    update_positions!(p, pmover.meshx, dt/2)
    update_velocities!(p, pmover, dt)
    update_positions!(p, pmover.meshx, dt/2)

    # Returns the square of the electric energy, computed in three different ways.
    return sum(p.wei .* pmover.phi), sum(pmover.∂phi.^2) * pmover.meshx.stop / p.nbpart, compute_int_E(p, pmover)
end

function PIC_step!(p, meshx, meshv, dt, dphidx)
    #=
    PIC step when supplying the derivative of the electrostatic potential
    =# 
    # Use a 3-step splitting, to be of order 2.
    # We interpolate the derivative of the potential, defined on a grid, to obtain an approximation
    # of its value at each particle position.
    for ipart = 1:p.nbpart
        idxgridx = Int64(fld(p.x[ipart], meshx.step)) + 1
        t = (p.x[ipart] - (idxgridx-1) * meshx.step) / meshx.step
        p.v[ipart] -= dt/2 * (dphidx[idxgridx] * (1-t) + dphidx[idxgridx < meshx.len ? idxgridx+1 : 1] * t)
        # Periodic boundary conditions in velocity
        if p.v[ipart] > meshv.stop
            p.v[ipart] -= meshv.stop - meshv.start
        elseif p.v[ipart] < meshv.start
            p.v[ipart] += meshv.stop - meshv.start
        end
    end
    update_positions!(p, meshx, dt)
    for ipart = 1:p.nbpart
        idxgridx = Int64(fld(p.x[ipart], meshx.step)) + 1
        t = (p.x[ipart] - (idxgridx-1) * meshx.step) / meshx.step
        p.v[ipart] -= dt/2 * (dphidx[idxgridx] * (1-t) + dphidx[idxgridx < meshx.len ? idxgridx+1 : 1] * t)
        # Periodic boundary conditions in velocity
        if p.v[ipart] > meshv.stop
            p.v[ipart] -= meshv.stop - meshv.start
        elseif p.v[ipart] < meshv.start
            p.v[ipart] += meshv.stop - meshv.start
        end
    end



    return sqrt.(sum(dphidx.^2) * meshx.stop / p.nbpart)
end

"""
compute rho, charge density (ie int f dv)

Projection on mesh is done here.
"""
function compute_rho!(p, pmover)
    nx = pmover.meshx.len
    dx = pmover.meshx.step
    pmover.ρ .= 0.0
 
    @inbounds @simd for ipart=1:p.nbpart
        idxonmesh = Int64(fld(p.x[ipart], dx)) + 1
        t = (p.x[ipart] - (idxonmesh-1) * dx) / dx
        if idxonmesh > pmover.meshx.len
            idxonmesh -= pmover.meshx.len
        elseif idxonmesh < 1
            idxonmesh += pmover.meshx.len
        end
        pmover.ρ[idxonmesh] += p.wei[ipart] * (1-t)
        pmover.ρ[idxonmesh < nx ? idxonmesh+1 : 1] += p.wei[ipart] * t   
    end
    pmover.ρ ./= dx
 
    pmover.ρ .-= sum(pmover.ρ) * dx / pmover.meshx.stop
end

function compute_phi!(pmover)
    dx = pmover.meshx.step
    L = pmover.meshx.stop
    
    pmover.phi_grid .= pmover.poisson_matrix \ pmover.ρ
    pmover.phi_grid .-= sum(pmover.phi_grid) * dx / L
end

function compute_E(pmover)
    return -(circshift(pmover.phi_grid, 1) .- circshift(pmover.phi_grid, -1)) ./ (2*pmover.meshx.step)
end

function compute_int_E(p, pmover)
    compute_rho!(p, pmover)
    compute_phi!(pmover)
    E = compute_E(pmover)
    return sum(E.^2) * pmover.meshx.step
end