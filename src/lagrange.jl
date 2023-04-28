import SemiLagrangian: Lagrange, interpolate!, get_order

export BSLLagrange

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct BSLLagrange <: AbstractMethod

    order:: Int
    interp:: Lagrange
    
    BSLLagrange( order ) = new( order, Lagrange(order, Float64))

end

"""
$(SIGNATURES)
"""
function advection!(fp, fi, mesh, interp::Lagrange, v, dt)

    dec = - v * dt / mesh.step
    decint = floor(Int, dec)
    value = dec - decint
    
    if get_order(interp) % 2 == 0 && value > 0.5
        value -= 1
        decint += 1
    end
        
    precal = [fct(value) for fct in interp.tabfct]    
    
    interpolate!(fp, fi, decint, precal, interp)
    
end
