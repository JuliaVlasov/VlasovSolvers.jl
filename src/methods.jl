abstract type AbstractMethod end

export BSL

struct BSL <: AbstractMethod

    p :: Int

end

export Fourier

struct Fourier <: AbstractMethod
    
    kx :: Vector{Float64}
    kv :: Vector{Float64}
      
end

