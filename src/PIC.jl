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
    x #list of the positions
    v #list of the velocities
    wei #list of the weights of the particules
    nbpart #nulber of particules
end


function samples(nsamples, kx, α::Float64, μ::Float64, β::Float64)
    #=
    Sample distribution defined by 1+\alpha \cos(kx * x) in space and a gaussian of variance beta and mean μ in velocity.
    =#
    x0 = Array{Float64}(undef, nsamples)
    seq = SobolSeq(1)
    for i in 1:nsamples
        r = Sobol.next!(seq)
        resx = find_zero(x-> x + α/kx * sin(kx*x) - r[1]*2π/kx, 0.0)
        x0[i] = resx
    end
    v0 = rand(Normal(μ, 1/√β), nsamples)
    wei = 2π/kx/nsamples .* ones(Float64, nsamples)
    return x0, v0, wei
end

function update_positions!(p, mesh, dt)
    for i = eachindex(p.x)
        p.x[i] += p.v[i] * dt 
        # Periodic boundary conditions
        if p.x[i] > mesh.stop
            p.x[i] -= mesh.stop
        elseif p.x[i] < mesh.start
            p.x[i] += mesh.stop
        end
    end
end

"""
update particle velocities vp (phi_v)
"""
function update_velocities!(p, n, dt, kx)
    S, C = compute_S_C(p, n, kx)
    der_phi = similar(p.x)
    phi = similar(p.x)
    for posidx = eachindex(p.x)
        der_phi[posidx] = 0
        phi[posidx] = 0
        for k = eachindex(C)
            der_phi[posidx] += 1/(pi*k) * (-sin(k*kx*p.x[posidx]) * C[k] + cos(k*kx*p.x[posidx]) * S[k])
            phi[posidx] += 1/(pi*k^2*kx) * (cos(k*kx*p.x[posidx]) * C[k] + sin(k*kx*p.x[posidx]) * S[k])
        end
    end
    p.v .-= dt .* der_phi
    return phi, der_phi
end

""" S[k] = \\sum_l=1^n {\\beta_l * sin(k kx x_l)} et C[k] = \\sum_l=1^n {\\beta_l * cos(k kx x_l)}
utile pour le calcul de la vitesse et du potentiel electrique phi"""
function compute_S_C(p, n, kx)
    S = Array{Float64}(undef, n)
    C = Array{Float64}(undef, n)
    for k = eachindex(S)
        S[k] = 0
        C[k] = 0
        for l = eachindex(p.x)
            S[k] += p.wei[l] * sin((k * kx) * p.x[l])
            C[k] += p.wei[l] * cos((k * kx) * p.x[l])
        end
    end
    return S, C
end

"""
compute rho, charge density (ie int f dv)
"""
function compute_rho(p, m)

    nx = m.len
    dx = m.step
    ρ = zeros(nx)
 
    for ipart=1:p.nbpart
        idxonmesh = Int64(fld(p.x[ipart], dx)) + 1
        t = (p.x[ipart] - (idxonmesh-1) * dx) / dx
        ρ[idxonmesh] += p.wei[ipart] * (1-t)
        ρ[idxonmesh < nx ? idxonmesh+1 : 1] += p.wei[ipart] * t   
    end
    ρ ./= dx
 
    ρ_tot  = sum(ρ) * dx / m.stop
    return ρ, ρ_tot
end

function compute_phi(rho, rho_total, mesh)
    dx = mesh.step
    L = mesh.stop
    nx = mesh.len
    matrix_poisson = spdiagm(  -1 => -ones(Float64,nx-1),
                                0 => 2*ones(Float64,nx),
                                1 => -ones(Float64,nx-1))
    matrix_poisson[1, nx] = -1
    matrix_poisson[nx, 1] = -1
    phi = matrix_poisson \ ((rho .- rho_total) .* dx^2 )
        # - rho_total car le monde est circulaire, sol périodique
    phi_total = sum(phi) * dx / L
    return phi .- phi_total # - phi_total densité car centré en 0
end

function compute_E(phi, mesh)
    return -(circshift(phi, 1) .- circshift(phi, -1)) ./ (2*mesh.step)
end

function compute_int_E(p, mesh)
    rho, rho_t = compute_rho(p, mesh)
    phi = compute_phi(rho, rho_t, mesh)
    E = compute_E(phi, mesh)
    s = 0.0
    for ide = eachindex(E)
        s += E[ide]^2
    end
    return s * mesh.step
end


function PIC_step(p, mesh, dt)
    energy_cine = 1 / 2 * sum(p.wei .* (p.v.^2))

    # Use a 3-step splitting, to be of order 2.
    update_velocities!(p, 1, dt/2, 2π/mesh.stop) # 1 because we are in HMF framework
    update_positions!(p, mesh, dt)
    phi_t, der_phi_t = update_velocities!(p, 1, dt/2, 2π/mesh.stop) # 1 because we are in HMF framework

    return sum(p.wei .* phi_t), compute_int_E(p, mesh), sum(der_phi_t.^2) * mesh.stop / length(p.x)
end