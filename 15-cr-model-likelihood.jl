# individual likelihood
function ill(θ::Vector{<:Number}, obs::Observation)
    d = obs.dose
    R = obs.response

    D = Multinomial(1, η(d, θ))
    
    return loglikelihood(D, R)
end


# negative log-likelihood function
function nll(ξ::Design, θ::Vector{<:Number}, obs::Vector{Observation})
    nllv = 0
    
    for i in eachindex(obs)
        nllv += ill(θ, obs[i])
    end

    return -nllv - 0.5 * log(det(FIM(ξ, θ) + 1e-10 .* I(4)))
end


# MLEs
function mle(ξ::Design, obs::Vector{<:Observation}, S::Number, N::Int64, L::Vector{<:Number}, U::Vector{<:Number})
    # defining negative log-likelihood function 
    f(θ) = nll(ξ, θ, obs)

    # finding MLEs by PSO
    swarm = Swarm(S, N, L, U, (0.9, 0.4), 1.25, 2.5, 0.5, 1e-16, f);
    θmle = optimize_mle!(f, swarm)

    return θmle
end