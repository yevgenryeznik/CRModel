# ======================================================================
#
# In this script, we define a continuation-ratio (CR) model used to
# model dose-efficacy-toxicity response in dose-finding studies.
#
# The response is modeled as a four-parameter vector-function: 
# η(x, θ) = (π0(x, θ), π1(x, θ), π2(x, θ))', where
#   - x is a given dose;
#   - θ = (θ[1], θ[2], θ[3], θ[4])' is a vector of model parameters.
#
# Here,  
#   - π0(x, θ) is a probability of "no effect";
#   - π1(x, θ) is a probability of "effect without toxicity";
#   - π2(x, θ) is a probability of "toxicity".
#
# ======================================================================


# sigmoid function used to evaluate CR model (values of η(x, θ) components) 
sigmoid(x::Number) = 1/(1 + exp(x))


# CR model evaluated at a single dose, given a set of parameters θ
function η(x::Number, θ::Vector{<:Number})
    # probability of "no effect"
    π0 = sigmoid(θ[1] + θ[2]*x)*sigmoid(θ[3] + θ[4]*x)

    # probability of "effect without toxicity"
    π1 = (1 - sigmoid(θ[1] + θ[2]*x))*sigmoid(θ[3] + θ[4]*x)

    # probability of "toxicity"
    π2 = 1 - sigmoid(θ[3] + θ[4]*x)

    return [π0, π1, π2]
end


# CR model evaluated at multiple doses, given a set of parameters θ
function η(x::Vector{<:Number}, θ::Vector{<:Number})
    return [η(x[i], θ) for i in eachindex(x)]
end


# a Jacobian of η(x, θ) with respect to θ, given a single dose x
function ∇η(x::Number, θ::Vector{<:Number})
    return jacobian(p -> η(x, p), θ)[1]
end


# a Jacobian of η(x, θ) with respect to θ, given multiple doses x
function ∇η(x::Vector{<:Number}, θ::Vector{<:Number})
    return [∇η(x[i], θ) for i in eachindex(x)]
end