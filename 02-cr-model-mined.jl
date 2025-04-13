# ======================================================================
#
# In this script, we define functions to calculate/estimate 
# a "Minimum Efficacious Dose" (MinED).
#
# On a continuous dose space [dmin, dmax], MinED is obtained by solving 
# an equation: π0(x, θ) = Δ, where Δ is a probability of "no efficacy".
#
# On a discrete dose space D, one may use  the following two definitions:
#   - rule #1: MinED is defined as a dose d* ∈ D, at which 
#     |π0(d*, θ) - Δ| has the lowest value among all d ∈ D.
#
#   - rule #2: MinED is defined as the lowest dose d* ∈ D, at which 
#     Δ ≥ π0(d*, θ).
#
# ======================================================================

# logit function, y = logit(x), defined as an inverse function to the 
# logistic function exp(x)/(1+exp(x)) =>
# y = log(x/(1-x)). The function used to calculate "Maximum Tolerated Dose" (MTD)
logit(x::Number) = log(x/(1-x))


# MinED on a continuous dose space [dmin, dmax] can be defined as
#   - a dose, at which probability of "neutrality" equals to Δ: rule #1 
#     (in this case, default Δ =0.2);
#   - a dose, at which probability of "efficcy without toxicity" equals to Δ: rule #2
#     (in this case, default Δ =0.6);
function mined(θ::Vector{<:Number}, Δ::Number = 0.2; rule::Int64 = 1)
    @assert rule ∈ [1, 2] "`rule` value must be in {1, 2, 3, 4}!"

    if rule == 1
        # MinED is a solution of π0(x, θ)} = Δ ⇔ g(x, θ) - Δ = 0
        g(x::Number, θ::Vector{<:Number}) = sigmoid(θ[1] + θ[2]*x)*sigmoid(θ[3] + θ[4]*x) - Δ

        # setting NonlinearProblem
        problem = NonlinearProblem(g, 0, θ)

        # solving
        solution = NonlinearSolve.solve(problem) 
        
        # extracting MinED
        MinED = solution.u

        # since MinED(θ) is a function of θ, one can calculate its gradient, ∇MinED(θ).
        # to evaluate ∇MinED(θ), we need to evaluate:
        #   - ∂g(MinED, θ)/∂x
        #   - ∂g(MinED, θ)/∂θ 
        # where g(MinED, θ) = 0
        ∂g∂x = Zygote.gradient(x -> g(x, θ), MinED)[1] 
        ∂g∂θ = Zygote.gradient(θ -> g(MinED, θ), θ)[1] 
        
        # then, ∇MinED(θ) is given by
        ∇MinED = (1/∂g∂x) .* ∂g∂θ

        # returning MED and ∇MinED
        return MinED, ∇MinED
    else
        # MinED is obtained by solving exp(θ[1] + θ[2])/(1 + exp(θ[1] + θ[2]x)) = Δ,
        # where Δ is a probability of efficacy without toxicity:
        # Δ ⇔ exp(θ[1] + θ[2]*x)/(1+exp(θ[1] + θ[2]*x)) = Γ ⇔ θ[1] + θ[2]*x = logit(Δ)
        
        # defining MTD, given Γ
        MinED = (logit(Δ) - θ[1])/θ[2]

        # calculating ∇mtd(θ)
        ∇MinED = [-1/θ[2], (θ[2] - logit(Δ))/θ[2]^2, 0, 0]

        # returning MTD
        return MinED, ∇MinED
    end
    

    
end


# MinED on a discrete dose D = {d[1], d[2], ..., d[n]}: rule #1
# default value of Δ is set to 0.2 (i.e., probability of "no effect" equals to 20%)
function mined(θ::Vector{<:Number}, D::Vector{<:Number}, Δ::Number = 0.2; rule::Int64 = 1)
    @assert rule ∈ [1, 2, 3, 4] "`rule` value must be in {1, 2, 3, 4}!"

    if rule in [1, 2]
        # responses corresponding to the "neutral" curve 
        y = [η(d, θ)[1] for d in D]
    else 
        # responses corresponding to the "efficacy without toxicyt" curve 
        y = [exp(θ[1] + θ[2] * d) / (1 + exp(θ[1] + θ[2] * d)) for d in D]
    end
    
    if rule in [1, 3]                      # MinED by rule #1 and #3
        _, i = findmin(abs.(y .- Δ))
        return D[i]
    elseif rule == 2                       # MinED by rule #2
        for d in D
            if Δ >= η(d, θ)[1]
                return d
            end
        end
    else                                   # MinED by rule #4
        for i in eachindex(y)
            if Δ <= y[i]
                return D[i]
            end
        end
    end
end





