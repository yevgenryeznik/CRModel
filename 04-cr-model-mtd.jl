# ======================================================================
#
# In this script, we define functions to calculate/estimate 
# a "Maximum Tolerated Dose" (MTD).
#
# On a continuous dose space [dmin, dmax], MTD is obtained by solving 
# an equation: π2(x, θ) = Γ, where Γ is a probability of "toxicity".
#
# On a discrete dose space D, one may use  the following two definitions:
#   - rule #1: MTD is defined as a dose d* ∈ D, at which 
#     |π2(d*, θ) - Γ| has the lowest value among all d ∈ D.
#
#   - rule #2: MTD is defined as the highest dose d* ∈ D, at which 
#     Γ ≤ π2(d*, θ).
#
# ======================================================================


# logit function, y = logit(x), defined as an inverse function to the 
# logistic function exp(x)/(1+exp(x)) =>
# y = log(x/(1-x)). The function used to calculate "Maximum Tolerated Dose" (MTD)
logit(x::Number) = log(x/(1-x))


# function calculates a "Maximum Tolerated Dose" (MTD)
# a search is conducted on an interval of doses
function mtd(θ::Vector{<:Number}, Γ::Number = 0.2)
    # MTD is obtained by solving π2(x, θ) = Γ, where Γ is a probability of toxicity:
    # π2(x, θ) = Γ ⇔ exp(θ[3] + θ[4]*x)/(1+exp(θ[3] + θ[4]*x)) = Γ ⇔ θ[3] + θ[4]*x = logit(Γ)
    
    # defining MTD, given Γ
    MTD = (logit(Γ) - θ[3])/θ[4]

    # calculating ∇mtd(θ)
    ∇MTD = [0, 0, -1/θ[4], (θ[3] - logit(Γ))/θ[4]^2]

    # returning MTD
    return MTD, ∇MTD
end


# function calculates a "Maximum Tolerated Dose" (MTD)
# a search is conducted in a set of discrete doses
function mtd(θ::Vector{<:Number}, D::Vector{<:Number}, Γ::Number = 0.2; rule::Int64 = 1)
    @assert rule ∈ [1, 2] "`rule` value must be in {1, 2}!"
    
    # responses corresponding to the "neutral" curve 
    y = [η(d, θ)[3] for d in D]
    
    
    if rule == 1                      # MTD by rule #1
        _, i = findmin(abs.(y .- Γ))
        return D[i]
    else                              # MTD by rule #2
        for d in sort(D, rev=true)
            if Γ >= η(d, θ)[3]
                return d
            end
        end
    end
end