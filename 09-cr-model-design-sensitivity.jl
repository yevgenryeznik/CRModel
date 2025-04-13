# ======================================================================
#
# In this script, we define functionality to deal with sensitivity functions
# used to test designs' optimality.
#
# A sensitivity function is defined as a directional derivative of 
# an optimality criterion Ψ at ξ in the direction of x (a single-point design)
#
# ======================================================================


# sensitivity function of ΨD
function ψD(x::Number, ξ::Design, θ::Vector{<:Number})
    # number of parameters
    m = length(θ)
    
    # Π = η(x, θ) = [π1(x, θ), π2(x, θ), π3(x, θ)]
    Π = η(x, θ)

    # elemental Information Matrix (EIM) for a multinomial distribution:
    # Y ~ Multinomial(1, Π)
    ν = diagm(1 ./ Π) - ones(3, 3)

    # Fisher information matrix (FIM)
    M = FIM(ξ, θ)

    # inverse FIM(ξ, θ)
    D = (M + 1e-10 .* I(m))\I(m)

    return tr(ν*∇η(x, θ)*D*∇η(x, θ)') - m    
end


# sensitivity function of Ψbod
function ψobd(x::Number, ξ::Design, θ::Vector{<:Number})
    # number of parameters
    m = length(θ)
    
    # Unit Information Matrix (UIM)
    Ix = μ(x, θ)

    # Fisher information matrix (FIM)
    M = FIM(ξ, θ)

    # inverse FIM(ξ, θ)
    D = (M + 1e-10 .* I(m))\I(m)
    
    # ∇OBD(θ)
    _, ∇ = obd(θ)
    
    # ∇OBD * ∇OBDᵀ
    ∇∇T = ∇ * transpose(∇)
    
    return tr(Ix*D*∇∇T*D) - tr(∇∇T*D)    
end


# sensitivity function of Ψmtd
function ψmtd(x::Number, ξ::Design, θ::Vector{<:Number}, Γ::Number = 0.2)
    # number of parameters
    m = length(θ)
    
    # Unit Information Matrix (UIM)
    Ix = μ(x, θ)

    # Fisher information matrix (FIM)
    M = FIM(ξ, θ)

    # inverse FIM(ξ, θ)
    D = (M + 1e-10 .* I(m))\I(m)
    
    # ∇MTD(θ)
    _, ∇ = mtd(θ, Γ)
    
    # ∇MTD * ∇MTDᵀ
    ∇∇T = ∇ * transpose(∇)
    
    return tr(Ix*D*∇∇T*D) - tr(∇∇T*D)    
end