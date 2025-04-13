# ======================================================================
#
# In this script, we define functionality to calculate designs efficiencies
# with respect to D- and c-optimal designs.
#
# ======================================================================


# function calculates D-efficiency
function deff(ξ::Design, ξD::Design, θ::Vector{<:Number})
    M = FIM(ξ, θ)
    MD = FIM(ξD, θ)

    return round((det(M)/det(MD))^0.25, digits = 2)
end


# function calculates c-efficiency (OBD)
function ceff(ξ::Design, ξc::Design, θ::Vector{<:Number})
    return round(Ψobd(ξc, θ)/Ψobd(ξ, θ), digits = 2)
end