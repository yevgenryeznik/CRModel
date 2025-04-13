# ==============================================================================
#
# In this script, we define a deisgn based on Random Walk Rule
# following to Ivanova A (2003)
#
# ============================================================================

function random_walk(nsbj::Int64, D::Vector{<:Number},θ::Vector{<:Number})
    i = 1
    R = zeros(Int64, 3, nsbj)
    d = zeros(Float64, nsbj)
    nD = length(D)

    for j in 1:nsbj
        d[j] = D[i]
        p = η(d[j], θ)
        R[:, j] = rand(Multinomial(1, p))
        if R[1, j] == 1
            i = i == nD ? i : i + 1
        elseif R[3, j] == 1
            i = i == 1  ? i : i - 1
        else
            i += 0
        end
    end

    return d
end 


function random_walk(D::Vector{<:Number}, θ::Vector{<:Number})
    # number of doses in a discrete dose set
    l = length(D)

    # auxilliary sequence
    λ = zeros(Float64, l)
    λ[2:l] = map(2:l) do j
        η(D[j-1], θ)[1]/η(D[j], θ)[3]    
    end
    λ[1] = 1 / (1 + sum([prod(λ[2:k]) for k in 2:l]))

    # steady-state allocation probabilities
    Π = [prod(λ[1:i]) for i in 1:l]
    Π = Π ./ sum(Π)
    
    return Design(D, Π)
end