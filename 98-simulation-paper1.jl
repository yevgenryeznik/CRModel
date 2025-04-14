# here we simulate results for the paper
#   "Nature-inspired Metaheuristics for Finding 
#   Efficient and Flexible Clinical Trials."

using CairoMakie
using DataFrames
using Distributions
using FLoops
using Latexify
using LaTeXStrings
using LinearAlgebra
using Match
using NonlinearSolve
using Pipe
using Printf
using Random: seed!
using Zygote


# including necessary scripts
include("01-cr-model-definition.jl")
include("02-cr-model-mined.jl")
include("03-cr-model-obd.jl")
include("04-cr-model-mtd.jl")
include("05-cr-model-visualize.jl")
include("06-cr-model-design.jl")
include("07-cr-model-design-fim.jl")
include("08-cr-model-design-optimality.jl")
include("09-cr-model-design-sensitivity.jl")
include("10-cr-model-design-pso.jl")
include("11-cr-model-design-check-optimality.jl")
include("12-cr-model-design-eff.jl")
include("13-cr-model-design-rwr.jl")
include("14-cr-model-observation.jl")
include("15-cr-model-likelihood.jl")

θ = [-3.5, 1, -6, 0.72]

# doses from the Alam paper
D = collect(0:2.5:10)
fig, _ = plot_cr_model(θ, 0, 10.0, add_bounds = false)
fig

MinED = round(mined(θ)[1], digits = 2)
OBD = round(obd(θ)[1],  digits = 2)
MTD = round(mtd(θ)[1], digits = 2)

# calculating optimal designs

D = 20
# lower bounds
L = [
    [  zeros(D ÷ 2); -ones(D ÷ 2)],
    [  zeros(D ÷ 2); -ones(D ÷ 2)],
    [  zeros(D ÷ 2); -ones(D ÷ 2)],
    [  zeros(D ÷ 2); -ones(D ÷ 2)]
]

# upper bounds
U = [
    [ 10*ones(D ÷ 2);  ones(D ÷ 2)],
    [MTD*ones(D ÷ 2);  ones(D ÷ 2)],
    [ 10*ones(D ÷ 2);  ones(D ÷ 2)],
    [MTD*ones(D ÷ 2);  ones(D ÷ 2)]
]

# objective functions
fD(P::Vector{<:Number}) = ofvD(P, θ)
fOBD(P::Vector{<:Number}) = ofvOBD(P, θ)

ξ = Design[]
for d in 1:4
    ofv = d in [1, 2] ? fD : fOBD
    seed!(3141592)
    swarm = Swarm(25, 700, L[d], U[d], (0.9, 0.4), 1.5, 2.5, 0.5, 1e-12, ofv)
    ξopt = optimize!(ofv, swarm)

    push!(ξ, ξopt)
end

# design I
ξ[1]

# design II
ξ[2]

# design III
ξ[3]

# design IV
ξ[4]

# checking optimality via sensitivity function (ψ)
d = 1 # design I 
# d = 2 # design II
# d = 3 # design III, 
# d = 4 # design IV

xlimits = [(0, 10), (0, MTD), (0, 10), (0, MTD)]
ylimits = [(-3, 0.5), (-2, 0.5), (-6, 1), (-8, 1)]
fig, ax = check_optimality(
    ξ[d], 
    θ, 
    d in [1, 2] ? x -> ψD(x, ξ[d], θ) : x -> ψobd(x, ξ[d], θ), 
    [xlimits[d][1], xlimits[d][2]],
    MinED,
    OBD,
    MTD
)
xlims!(ax, xlimits[d][1], xlimits[d][2]+0.5)
ylims!(ax, ylimits[d])
fig 


# simulating design for estimating operating characteristics
function fixed_optimal_design(ξ::Design, nsbj::Int64, θ::Vector{<:Number}, Lθ::Vector{<:Number}, Uθ::Vector{<:Number}; N = 50, max_iter = 1000, seed = 314159)    
    # simulating one adaptive design
    seed!(seed);
    
    # generating observations
    obs = generate_obs(nsbj, ξ, θ)

    # defining negative log-likelihood function 
    f(θ) = nll(θ, obs)

    # finding MLEs by PSO
    swarm = Swarm(N, max_iter, Lθ, Uθ, (0.9, 0.4), 1.5, 2.5, 0.5, 1e-12, f);
    θmle = optimize_mle!(f, swarm)

    return θmle, obs
end

# sample size
nsbj = 40

# number of simulations
nsim = 1000

ϵ = 1e-10
b = 20
Lθ = [-b, ϵ, -b, ϵ]
Uθ = [ b, b, -ϵ, b]

mle_θ = Vector[]
obs_ξ = Vector[]


for d in eachindex(ξ)
    push!(mle_θ, Vector[])
    push!(obs_ξ, Vector[])

    @time for s in 1:nsim
        mle_, obs_ = fixed_optimal_design(ξ[d], nsbj, θ, Lθ, Uθ, seed = 314159+2*s);
        push!(mle_θ[d], mle_)
        push!(obs_ξ[d], obs_)

        println("simulation $s is complete!")
    end    
end



function summarize(mle_θ, design)
    OBDhat = [obd(item)[1] for item in mle_θ]
    MTDhat = [mtd(item)[1] for item in mle_θ]

    keep = (0 .< MTDhat .< 10) .& (0 .< OBDhat .< 10)

    OBD = obd(θ)[1]    
    OBDbias = OBD - mean(OBDhat[keep])
    OBDsd = std(OBDhat[keep])
    OBDrmse = sqrt(OBDbias^2 + OBDsd^2)

    MTD = mtd(θ)[1]
    MTDbias = MTD - mean(MTDhat[keep])
    MTDsd = std(MTDhat[keep])
    MTDrmse = sqrt(MTDbias^2 + MTDsd^2)

    Pr_NoMLE = round(1-mean(keep), digits = 4)

    @pipe DataFrame(
        "MTD-Bias" => round(MTDbias, digits = 4),
        "MTD-SD"   => round(MTDsd, digits = 4),
        "MTD-RMSE" => round(MTDrmse, digits = 4),
        "OBD-Bias" => round(OBDbias, digits = 4),
        "OBD-SD"   => round(OBDsd, digits = 4), 
        "OBD-RMSE" => round(OBDrmse, digits = 4),
        "Pr(No MLE)" => Pr_NoMLE
    ) |>
    stack(_, 
        names(_), 
        variable_name = "characteristic", 
        value_name = design
    )
end

# operating characteristics of the 4 designs
innerjoin(
    summarize(mle_θ[1], "I"),
    summarize(mle_θ[2], "II"),
    summarize(mle_θ[3], "III"),
    summarize(mle_θ[4], "IV"),
    on = :characteristic
) 


