export Fourier

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct Fourier <: AbstractMethod
    
    kx   :: Vector{Float64}
    kv   :: Vector{Float64}

    function Fourier( meshx, meshv )
    
        nx  = meshx.len
        dx  = meshx.step
        Lx  = meshx.stop - meshx.start
        kx  = 2π/Lx .* [0:nx÷2-1;-nx÷2:-1]

        nv  = meshv.len
        dv  = meshv.step
        Lv  = meshv.stop - meshv.start
        kv  = 2π/Lv .* [0:nv÷2-1;-nv÷2:-1]

        new( kx, kv)

    end
  
end


"""
$(SIGNATURES)
"""
function advection_v!(fᵗ  :: Array{ComplexF64,2}, 
                      adv :: Fourier,
		      e   :: Vector{ComplexF64}, 
		      dt  :: Float64 )
    fft!(fᵗ, 1)
    @. fᵗ *= exp(-1im * dt * adv.kv * $transpose(e))
    ifft!(fᵗ, 1)

end

"""
$(SIGNATURES)
"""
function advection_x!( f   :: Array{ComplexF64,2}, 
                       adv :: Fourier,
		       e   :: Vector{ComplexF64}, 
		       v   :: Vector{Float64}, 
		       dt  :: Float64 )
    
    ev = exp.(-1im*dt * adv.kx * transpose(v))    
    
    fft!(f,1)
    f .= f .* ev
    dv = v[2]-v[1]
    ρ = dv * vec(sum(f,dims=2))  
    for i in 2:length(e)
        e[i] = -1im * ρ[i] ./ adv.kx[i]
    end
    e[1] = 0.0
    ifft!(f,1)
    ifft!(e)
    e .= real(e)
end


