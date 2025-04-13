# ======================================================================
#
# In this script, we define functions to calculate/estimate 
# an "Optimal Biological Dose" (OBD).
#
# On a continuous dose space [dmin, dmax], OBD is obtained by finding 
#     x* s.t. π1(x*, θ) = max{π1(x, θ)} for all x ∈ [dmin, dmax].
#
# On a discrete dose space D, OBD is obtained by finding 
#     d* s.t. π1(d*, θ) = max{π1(d, θ)} for all d ∈ D.
#
# ======================================================================

# function calculates a "Optimal Biological Optimal Dose" (OBD)
# a search is conducted on an interval of doses
function obd(θ::Vector{<:Number})
    # BOD is obtained by solving ∂/∂x{π1(x, θ)} = 0.
    # after simplification, the equation ∂/∂x{π1(x, θ)} = 0 is equivalent to g(x, θ) = 0, where:
    g(x::Number, θ::Vector{<:Number}) = θ[2]*(1 + exp(-θ[3] - θ[4]*x)) - θ[4]*(1 + exp(θ[1] + θ[2]*x))

    # setting NonlinearProblem
    problem = NonlinearProblem(g, 0, θ)

    # solving
    solution = NonlinearSolve.solve(problem) 
    
    # extracting OBD
    OBD = solution.u

    # since OBD(θ) is a function of θ, one can calculate its gradient, ∇OBD(θ).
    # to evaluate ∇OBD(θ), we need to evaluate:
    #   - ∂g(OBD, θ)/∂x
    #   - ∂g(OBD, θ)/∂θ 
    # where g(OBD, θ) = 0
    ∂g∂x = Zygote.gradient(x -> g(x, θ), OBD)[1] 
    ∂g∂θ = Zygote.gradient(θ -> g(OBD, θ), θ)[1] 
    
    # then, ∇bod(θ) is given by
    ∇OBD = (1/∂g∂x) .* ∂g∂θ

    # returning OBD
    return OBD, ∇OBD
end


# function calculates a "Optimal Biological Optimal Dose" (OBD)
# a search is conducted in a set of discrete doses
function obd(θ::Vector{<:Number}, D::Vector{<:Number})
    # responses corresponding to the "efficacy without toxicity" curve 
    y = [η(d, θ)[2] for d in D]
    _, i = findmax(y)

    return D[i]
end