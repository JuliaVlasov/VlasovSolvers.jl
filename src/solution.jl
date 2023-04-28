using RecipesBase

"""
$(TYPEDEF)

Data structure that stores the solution of the Vlasov problem of one dimension
in both physical space and phase space.

$(TYPEDFIELDS)
"""
struct VlasovSolution1D1V

    values :: Vector{Array{Float64, 2}}
    times :: Vector{Float64}
    energy :: Vector{Float64}
    saveat :: Float64

    function VlasovSolution1D1V( )

        values = Array{Float64, 2}[]
        times = Float64[]
        energy = Float64[]
        new( values, times, energy, 0.0)

    end

    function VlasovSolution1D1V( saveat :: Float64 )

        values = Array{Float64, 2}[]
        energy = Float64[]
        times = Float64[]
        new( values, times, energy, saveat)

    end
    
end


@recipe function plot(sol::VlasovSolution1D1V)

    title --> "Vlasov problem 1D1V"
    xlabel --> "time"
    ylabel --> "energy"
    
    sol.times, sol.energy
    

end
