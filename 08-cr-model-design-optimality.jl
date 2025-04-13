# ======================================================================
#
# In this script, we define functionality to calculate optimality criteria
# to find "optimal designs".
#
# ======================================================================


# D-optimal criteria for estimating model parameters
function ΨD(ξ::Design, θ::Vector{<:Number})
    # number of parameters
    m = length(θ)
    
    # FIM 
    M = FIM(ξ, θ)

    # returning D-optimality criterion
    return -log(det(M + 1e-10 .* I(m)))
end


# c-optimal criteria for estimating OBD
function Ψobd(ξ::Design, θ::Vector{<:Number})    
    # number of parameters
    m = length(θ)
    
    # FIM 
    M = FIM(ξ, θ)

    # inverse FIM(ξ, θ)
    D = (M + 1e-10 .* I(m))\I(m)

    # evaluating gradient of ∇OBD(θ)
    _, ∇OBD = obd(θ)

    return transpose(∇OBD) * D * ∇OBD
end


# c-optimal criteria for estimating MTD
function Ψmtd(ξ::Design, θ::Vector{<:Number}, Γ::Number = 0.2)    
    # number of parameters
    m = length(θ)
    
    # FIM 
    M = FIM(ξ, θ)

    # inverse FIM(ξ, θ)
    D = (M + 1e-10 .* I(m))\I(m)

    # ∇MTD(θ)
    _, ∇MTD = mtd(θ, Γ)

    return transpose(∇MTD) * D * ∇MTD
end