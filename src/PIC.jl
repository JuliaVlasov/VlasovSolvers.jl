using Sobol
using Roots

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
    Defines Runge-Kutta-Nystrom time integrator via its Butcher tableau,
    and holds some pre-allocated arrays used for the time integration only.
"""
struct rkn_order4
    a :: Array{Float64, 2}
    b̄ :: Vector{Float64}
    c :: Vector{Float64}
    b :: Vector{Float64}
    dt :: Float64
    fg ::       Array{Float64, 2}
    G ::        Array{Float64, 1}

    function rkn_order4(X, dt)
        # a, b̄, c, b correspond to the Butcher tableau of Runge-Kutta-Nystrom 3steps order4.
        a = [0.0        0.0       0.0; 
            (2-√3)/12   0.0           0.0; 
            0.0         √(3)/6      0.0]
        b̄ = [(5 - 3*√3)/24,     (3+√3)/12,  (1+√3)/24]
        c = [(3+√3)/6,          (3-√3)/6,   (3+√3)/6]
        b = [(3-2*√3)/12,       1/2,        (3+2*√3)/12]

        new(a .* dt^2, b̄ .* dt^2, c .* dt, b .* dt, dt, 
            zeros(Float64, length(X), 3),  # fg
            similar(X) # G
        )
    end
end


"""
Holds pre-allocated arrays
"""
struct ParticleMover
    Φ ::Vector{Float64}
    ∂Φ :: Vector{Float64}
    meshx
    kx
    K
    C ::Vector{Float64}
    S ::Vector{Float64}
    tmpcosk ::Vector{Float64}
    tmpsink ::Vector{Float64}
    tmpcoskimplicit ::Vector{Float64}
    tmpsinkimplicit ::Vector{Float64}
    poisson_matrix
    ρ ::Vector{Float64}
    Φ_grid ::Vector{Float64}
    idxonmesh :: Vector{Int64}
    idxonmeshp1 :: Vector{Int64}
    rkn :: rkn_order4
    dt :: Float64
    
    function ParticleMover(particles::Particles, meshx, K, dt; kx=1)
        Φ = similar(particles.x)
        ∂Φ = similar(particles.x)
        tmpcosk = similar(particles.x)
        tmpsink = similar(particles.x)
        nx = meshx.len
        np = particles.nbpart

        matrix_poisson = spdiagm(  -1 => .- ones(Float64,nx-1),
                                    0 => 2 .* ones(Float64,nx),
                                    1 => .- ones(Float64,nx-1))
        matrix_poisson[1, nx] = -1
        matrix_poisson[nx, 1] = -1

        matrix_poisson ./= meshx.step^2


        new(Φ, ∂Φ, meshx, kx, K, Vector{Float64}(undef, K), Vector{Float64}(undef, K), tmpcosk, tmpsink, similar(tmpcosk), similar(tmpsink), matrix_poisson, Vector{Float64}(undef, nx), Vector{Float64}(undef, nx), Vector{Int64}(undef, np), Vector{Int64}(undef, np), rkn_order4(particles.x, dt), dt)
    end
end



#==== Time steppers ====#
"""symplectic_RKN_order4!(X, V, F, rkn, kx)
    
    Advect (X, V) on a time step dt using symplectic Runge-Kutta-Nystrom method of order4 [Feng, Qin (2010), sect.7.3, p.327, scheme1].

    The equation satisfied by X is
    ```math
    \\frac{d^2 X(t)}{dt^2} = C(t)\\cos(X(t)) - S(t)\\sin(X(t))
    ```

    RKN method considers Ẋ = V as a variable, and updates both X and V.

    Args:
    - X: matrix of positions at time t_n
    - V: matrix of velocities at time t_n
    - F: values of initial condition at time t_0
    - rkn: rkn_order_4 struct, storing butcher tableau and pre-allocated arrays (holds the value of dt)
    - kx: 2π/L

    Updates X, V in place, and returns coefficients C, S at current time.

"""
function symplectic_RKN_order4!(p, pmover, kernel)
    @views begin
        for s=1:3
            pmover.rkn.G .= p.x .+ p.v .* pmover.rkn.c[s] .+ pmover.rkn.a[s, 1] .* pmover.rkn.fg[:, 1] .+ 
                                                             pmover.rkn.a[s, 2] .* pmover.rkn.fg[:, 2] .+ 
                                                             pmover.rkn.a[s, 3] .* pmover.rkn.fg[:, 3]
     
            kernel(pmover.rkn.fg[:, s], pmover.rkn.G, p, pmover; coeffdt=0)
        end
    
        p.x .+= pmover.dt .* p.v .+ pmover.rkn.b̄[1] .* pmover.rkn.fg[:, 1] .+ 
                                    pmover.rkn.b̄[2] .* pmover.rkn.fg[:, 2] .+ 
                                    pmover.rkn.b̄[3] .* pmover.rkn.fg[:, 3]
        p.v .+= pmover.rkn.b[1] .* pmover.rkn.fg[:, 1] .+ 
                pmover.rkn.b[2] .* pmover.rkn.fg[:, 2] .+ 
                pmover.rkn.b[3] .* pmover.rkn.fg[:, 3]
        
        kernel(pmover.∂Φ, p.x, p, pmover)
    end
end


"""strang_splitting!(particles, pmover, kernel)  
    
    Other method for advecting (X, V) on a time step. 

    Uses Verlet scheme (of order 2).

    Args:
    - particles: Particle struct
    - pmover: ParticleMover struct
    - kernel: how to compute the flow field.

    Updates X, V in place, and returns coefficients C, S at current time. 
"""
function strang_splitting!(particles, pmover, kernel)  
    pmover.Φ .= 0
    pmover.∂Φ .= 0
    particles.x .+= particles.v .* pmover.dt/2
    kernel(pmover.∂Φ, particles.x, particles, pmover)
    particles.v .+= pmover.dt .* pmover.∂Φ
    particles.x .+= particles.v .* pmover.dt/2
    kernel(pmover.∂Φ, particles.x, particles, pmover)
end

function strang_splitting_implicit!(particles, pmover, kernel)  
    pmover.Φ .= 0
    pmover.∂Φ .= 0
    κ = 1/2
    kernel(pmover.∂Φ, particles.x, particles, pmover; coeffdt=κ)
    particles.v .+= pmover.dt./2 .* pmover.∂Φ
    particles.x .+= particles.v .* pmover.dt
    kernel(pmover.∂Φ, particles.x, particles, pmover; coeffdt=κ)
    particles.v .+= pmover.dt./2 .* pmover.∂Φ
    kernel(pmover.∂Φ, particles.x, particles, pmover)
end


# ===== Kernel computations ==== #
"""
Compute -∂_x Φ[f](`x` + `pmover.dt*coeffdt*p.v`) and stores it in `dst`. Also updates `pmover.Φ`
"""
function kernel_poisson!(dst, x, p, pmover; coeffdt=0)
    dst .= 0
    
    if coeffdt > 0
        pmover.tmpcoskimplicit .= similar(pmover.tmpcosk)
        pmover.tmpsinkimplicit .= similar(pmover.tmpsink)
    end
    
    pmover.Φ .= 0

    for k=1:pmover.K        
        pmover.tmpcosk .= cos.(x .* (k * pmover.kx))
        pmover.tmpsink .= sin.(x .* (k * pmover.kx))

        if coeffdt > 0
            pmover.tmpcoskimplicit .= cos.((x .+ pmover.dt.*coeffdt.*p.v) .* (k*pmover.kx))
            pmover.tmpsinkimplicit .= sin.((x .+ pmover.dt.*coeffdt.*p.v) .* (k*pmover.kx))
        end
        
        pmover.C[k] = sum(pmover.tmpcosk .* p.wei)
        pmover.S[k] = sum(pmover.tmpsink .* p.wei)
        
        if coeffdt > 0 
            pmover.C[k] += -sum(p.v .* pmover.tmpsink .* p.wei) * k * 2π / pmover.meshx.stop * pmover.dt * coeffdt
            pmover.S[k] +=  sum(p.v .* pmover.tmpcosk .* p.wei) * k * 2π / pmover.meshx.stop * pmover.dt * coeffdt
        end

        if coeffdt > 0
            pmover.tmpcosk .= pmover.tmpcoskimplicit
            pmover.tmpsink .= pmover.tmpsinkimplicit
        end
        
        pmover.Φ .+= (pmover.C[k] .* pmover.tmpcosk .+ pmover.S[k] .* pmover.tmpsink) ./ k^2
        # The line below computes -∂Φ[f](`x`) and stores it to `dst`. 
        # Change the ".+=" to a ".-=" to have an attractive potential.
        dst .+= (pmover.C[k] .* pmover.tmpsink .- pmover.S[k] .* pmover.tmpcosk) ./ k 
    end
    
    dst ./= π
    pmover.Φ .*= pmover.meshx.stop / (2*π^2)
end

function kernel_gyrokinetic!(dst, x, p, pmover)
    dst .= .-x
end



# ===== Some quantities we can compute at each step ==== #
"""compute_electricalenergy²(p, pmover)

    Returns the square of the electrical energy
"""
function compute_electricalenergy²(p, pmover)
    return sum(p.wei .* pmover.Φ)
end


function compute_momentum(particles)
    return sum(particles.wei .* particles.v)
end

function compute_totalenergy²(particles, Eelec²)
    return 1/2 * (Eelec² .+ sum(particles.v.^2 .* particles.wei))
end



""" periodic_boundary_conditions!(p, pmover)

    Impose periodic boundary conditions in space.
"""
function periodic_boundary_conditions!(p, pmover)
    p.x[findall(>(pmover.meshx.stop), p.x)] .-= pmover.meshx.stop
    p.x[findall(<(0)                , p.x)] .+= pmover.meshx.stop
end



function PIC_step!(p::Particles, pmover::ParticleMover; kernel=kernel_poisson!)
    symplectic_RKN_order4!(p, pmover, kernel)
    # strang_splitting!(p, pmover, kernel)
    # strang_splitting_implicit!(p, pmover, kernel)
    
    periodic_boundary_conditions!(p, pmover)
    
    E² = compute_electricalenergy²(p, pmover)
    return E², compute_momentum(p), compute_totalenergy²(p, E²)
end


function PIC_step!(p, meshx, meshv, dt, dΦdx)
    #=
    PIC step when supplying the derivative of the electrostatic potential
    =# 
    # Use a 3-step splitting, to be of order 2.
    # We interpolate the derivative of the potential, defined on a grid, to obtain an approximation
    # of its value at each particle position.

    idxonmesh = Int64.(fld.(p.x, meshx.step)) .+ 1
    t = (p.x .- (idxonmesh.-1) .* meshx.step) ./ meshx.step
    idxonmesh[findall(i -> i > meshx.len,  idxonmesh)] .-= meshx.len
    idxonmesh[findall(i -> i < 1               ,  idxonmesh)] .+= meshx.len
    idxonmeshp1 = idxonmesh .+ 1
    idxonmeshp1[findall(i -> i > meshx.len,  idxonmeshp1)] .-= meshx.len
    p.v .-= dt ./ 2 .* (dΦdx[idxonmesh] .* (1 .- t) .+ dΦdx[idxonmeshp1] .* t)

    # for ipart = 1:p.nbpart
    #     idxgridx = Int64(fld(p.x[ipart], meshx.step)) + 1
    #     t = (p.x[ipart] - (idxgridx-1) * meshx.step) / meshx.step
    #     p.v[ipart] -= dt/2 * (dΦdx[idxgridx] * (1-t) + dΦdx[idxgridx < meshx.len ? idxgridx+1 : 1] * t)
    #     # # Periodic boundary conditions in velocity
    #     # if p.v[ipart] > meshv.stop
    #     #     p.v[ipart] -= meshv.stop - meshv.start
    #     # elseif p.v[ipart] < meshv.start
    #     #     p.v[ipart] += meshv.stop - meshv.start
    #     # end
    # end
    update_positions!(p, meshx, dt)
    # for ipart = 1:p.nbpart
    #     idxgridx = Int64(fld(p.x[ipart], meshx.step)) + 1
    #     t = (p.x[ipart] - (idxgridx-1) * meshx.step) / meshx.step
    #     p.v[ipart] -= dt/2 * (dΦdx[idxgridx] * (1-t) + dΦdx[idxgridx < meshx.len ? idxgridx+1 : 1] * t)
    #     # Periodic boundary conditions in velocity
    #     # if p.v[ipart] > meshv.stop
    #     #     p.v[ipart] -= meshv.stop - meshv.start
    #     # elseif p.v[ipart] < meshv.start
    #     #     p.v[ipart] += meshv.stop - meshv.start
    #     # end
    # end

    idxonmesh .= Int64.(fld.(p.x, meshx.step)) .+ 1
    t = (p.x .- (idxonmesh.-1) .* meshx.step) ./ meshx.step
    idxonmesh[findall(i -> i > meshx.len,  idxonmesh)] .-= meshx.len
    idxonmesh[findall(i -> i < 1               ,  idxonmesh)] .+= meshx.len
    idxonmeshp1 .= idxonmesh .+ 1
    idxonmeshp1[findall(i -> i > meshx.len,  idxonmeshp1)] .-= meshx.len
    p.v .-= dt ./ 2 .* (dΦdx[idxonmesh] .* (1 .- t) .+ dΦdx[idxonmeshp1] .* t)



    return sqrt.(sum(dΦdx.^2) * meshx.stop / p.nbpart)
end






















#= ================= #
    CLASSICAL PIC 
#  ================ =#

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
    p.x .+= p.v .* dt
    p.x[findall(x -> x >= mesh.stop,  p.x)] .-= mesh.stop - mesh.start
    p.x[findall(x -> x <  mesh.start, p.x)] .+= mesh.stop - mesh.start
end

"""
update particle velocities vp (Φ_v)
"""
function update_velocities!(p, pmover, dt)
    compute_S_C!(p, pmover)
    pmover.Φ .= 0
    pmover.∂Φ .= 0

    for k = 1:pmover.K
        pmover.tmpcosk .= cos.(k .* pmover.kx .* p.x)
        pmover.tmpsink .= sin.(k .* pmover.kx .* p.x)
        pmover.Φ   .+= (pmover.tmpcosk .* pmover.C[k] .+ pmover.tmpsink .* pmover.S[k]) ./ k^2
        pmover.∂Φ  .+= (.-pmover.tmpsink .* pmover.C[k] .+ pmover.tmpcosk .* pmover.S[k]) ./ k
    end
    pmover.Φ .*= pmover.meshx.stop / 2*π^2
    pmover.∂Φ ./= π
    p.v .-= dt .* pmover.∂Φ
end

""" S[k] = \\sum_l=1^n {\\beta_l * sin(k kx x_l)} et C[k] = \\sum_l=1^n {\\beta_l * cos(k kx x_l)}
utile pour le calcul de la vitesse et du potentiel electrique Φ"""
function compute_S_C!(p, pmover)
    pmover.S .= 0
    pmover.C .= 0
    
    for k = 1:pmover.K
        pmover.S[k] = sum(p.wei .* sin.(k .* 2π ./ pmover.meshx.stop .* p.x))
        pmover.C[k] = sum(p.wei .* cos.(k .* 2π ./ pmover.meshx.stop .* p.x))
    end
end





# ##### functions used to perform projection on grid. ##### #


"""
compute rho, charge density (ie int f dv)

Projection on mesh is done here.
"""
function compute_rho!(p, pmover)
    dx = pmover.meshx.step
    pmover.ρ .= 0.0
 
    pmover.idxonmesh .= Int64.(fld.(p.x, dx)) .+ 1
    t = (p.x .- (pmover.idxonmesh.-1) .* dx) ./ dx
    
    pmover.idxonmesh[findall(i -> i > pmover.meshx.len,  pmover.idxonmesh)] .-= pmover.meshx.len
    pmover.idxonmesh[findall(i -> i < 1               ,  pmover.idxonmesh)] .+= pmover.meshx.len
    
    pmover.idxonmeshp1 .= pmover.idxonmesh .+ 1
    pmover.idxonmeshp1[findall(i -> i > pmover.meshx.len,  pmover.idxonmeshp1)] .-= pmover.meshx.len
    
    pmover.ρ[pmover.idxonmesh]   .+= p.wei .* (1 .- t)
    pmover.ρ[pmover.idxonmeshp1] .+= p.wei .* t

    println((sum(pmover.ρ.^2), sum((p.wei .* t).^2), sum((p.wei .* (1 .-t)).^2)))

    pmover.ρ .-= sum(pmover.ρ) / pmover.meshx.stop
    return pmover.ρ
end


function compute_Φ!(pmover) # solving poisson equation with FD solver
    dx = pmover.meshx.step
    L = pmover.meshx.stop
    
    pmover.Φ_grid .= pmover.poisson_matrix \ pmover.ρ
    pmover.Φ_grid .-= sum(pmover.Φ_grid) * dx / L
end

function compute_int_E(p, pmover) # energy computed with a projection on the grid
    compute_rho!(p, pmover)
    compute_Φ!(pmover)
    E = compute_E(pmover)
    return sum(E.^2) * pmover.meshx.step
end

function compute_E(pmover)
    return -(circshift(pmover.Φ_grid, 1) .- circshift(pmover.Φ_grid, -1)) ./ (2*pmover.meshx.step)
end


# Carybe function
"""
compute rho, charge density (ie int f dv)
"""
function compute_rho( p, m )

   L = m.stop
   nx = m.len
   dx = m.step
   rho = zeros(nx)

   for ipart=1:p.nbpart
        xp = p.x[ipart]
        i = Int(floor(Real(xp) / dx) + 1) #number of the cell with the particule
        dum = p.wei[ipart] / dx
        mx_i = (i-1) * dx
        mx_j = i * dx
        a1 = (mx_j-xp) * dum
        a2 = (xp-mx_i) * dum
        if i < 1
            i += nx 
        elseif i > nx
            i -= nx
        end
        rho[i]   += a1
        rho[i < nx ? i+1 : i+1 - nx] += a2
   end

#    rho[nx] += rho[1]
#    rho[1] = rho[nx]

   rho ./= dx

   rho_total  = sum(rho) * dx / L
   return rho, rho_total
end