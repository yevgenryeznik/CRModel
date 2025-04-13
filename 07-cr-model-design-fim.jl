# ======================================================================
#
# In this script, we define functionality to calculate 
# Fisher Information Matrix (FIM). 
#
# ======================================================================


# Unit Information Matrix (UIM)
function μ(x::Number, θ::Vector{<:Number})
    # Π = η(x, θ) = [π0(x, θ), π1(x, θ), π2(x, θ)] 
    Π = η(x, θ) .+ 1e-10

    # elemental information matrix (EIM) for a multinomial distribution:
    # Y ~ Multinomial(1, Π)
    ν = diagm(1 ./ Π) - ones(3, 3)

    return ∇η(x, θ)' * ν * ∇η(x, θ)
end


# normalized Fisher Information Matrix for design `ξ`
function FIM(ξ::Design, θ::Vector{<:Number})
    # number of parameters
    m = length(θ)

    # FIM
    M = zeros(Float64, m, m)
    for (x, w) in zip(ξ.x, ξ.w)
        M += w .* μ(x, θ)
    end

    return M 
end