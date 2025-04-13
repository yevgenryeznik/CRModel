# defining a type, representing a single observation: a pair (dose, response)
struct Observation
    dose::Number
    response::Vector{Int64}   
end


function Base.show(io::IO, obs::Union{Observation, Vector{Observation}})
    if typeof(obs) == Observation
        println(io, "(dose: $(obs.dose), response: $(obs.response))")
    else
        for i in eachindex(obs)
            println(io, "(dose: $(obs[i].dose), response: $(obs[i].response))")
        end
    end

    return nothing
end


# function implements a Truncated Multinomial Randomization (TMD) procedure
function tmd(w::Vector{<:Number}, nsbj::Int64)
    target = Int.(floor.(w .* nsbj))
    ρ = target ./ sum(target)

    trt = zeros(Int64, nsbj)
    N = zeros(Int64, length(ρ))
    j = 1
    trtid = collect(eachindex(N))
    while any(N .< target)
        trt[j] = rand(trtid[N .< target])
        N[trt[j]] += 1
        j += 1
    end
    for i in j:nsbj
        trt[i] = rand(trtid)
    end
    return trt
end


# function generates response from the CR model, given a single dose `x`
function response(x::Number, θ::Vector{<:Number})
    D = Multinomial(1, η(x, θ))
    R = rand(D)

    return Observation(x, R)
end


# function generate observations from the CR model according to the design ξ
function generate_obs(nsbj::Int64, ξ::Design, θ::Vector{<:Number})
    x = ξ.x
    w = ξ.w
    
    dose = x[tmd(w, nsbj)]
    obs = [response(d, θ) for d in dose]

    return obs
end