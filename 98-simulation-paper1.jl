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


