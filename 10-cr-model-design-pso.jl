# ==============================================================================
#
# In this script, we define functionality to perform Particle Swarm Optimization
# to find optimal designs
#
# ============================================================================


# function to transform coordinates of a particle, 
# corresponding to weights, to values between 0 and 1
# and their sum equals to 1. 
normalize(u::Vector{<:Number}) = u.*u ./ sum(u.*u)


# a type representing a Particle
mutable struct Particle
    P::Vector{<:Number}     # particle position  
    V::Vector{<:Number}     # particle velocity
    B::Vector{<:Number}     # best particle position
    f::Float64              # OFV
    F::Float64              # best OFV
end


# a type representing a Swarm
mutable struct Swarm
    p::Vector{Particle}     # an array of particles
    S::Int64                # swarm size
    N::Int64                # number of iterations
    L::Vector{<:Number}     # lower bounds for design points
    U::Vector{<:Number}     # upper bounds for design points
    G::Vector{<:Number}     # best particle position in a swarm
    g::Float64              # best OFV in a swarm
    W::Vector{<:Number}     # inertia coefficients vs. iteration
    γ::Number               # relaxation parameter 
    C1::Vector{<:Number}    # cognitive parameters vs. iteration
    C2::Vector{<:Number}    # social parameters vs. iteration
    stop::Bool              # stopping flag
    ε::Float64              # tolerance

    # Swarm constructor
    function Swarm(
        S::Int64,
        N::Int64,
        L::Vector{<:Number}, 
        U::Vector{<:Number}, 
        w::Tuple{Number, Number},
        γ::Number, 
        c1::Number,
        c2::Number,
        ε::Float64,
        ofv::Function)

        # setting up the strategic parameters
        W  = [w[2] + ((N-i)/(N-1))^γ*(w[1] - w[2]) for i in 1:N] 
        C1 = [c1 - (i-1)/(N-1)*(c1-c2) for i in 1:N]
        C2 = [c2 + (i-1)/(N-1)*(c1-c2) for i in 1:N]

        p = Particle[]
        G = Number[]
        g = Inf
        
        for _ in 1:S
            # initializing particle position
            P = [low == upp ? low : rand(Uniform(low, upp)) for (low, upp) in zip(L, U)]

            # initializing velocities
            V = [rand(Uniform(-1, 1)) for _ in eachindex(P)]

            # evaluating objective function
            f = ofv(P)

            if (f < g)
                G = P
                g = f
            end
            particle = Particle(P, V, P, f, f)
        
            push!(p, particle)
        end

        return new(p, S, N, L, U, G, g, W, γ, C1, C2, false, ε)
    end
end


# functions used to convert a particle position into a designs
function position2design(P::Vector{Float64})
    # number of design points
    D = length(P) ÷ 2

    # design points
    x = P[1:D]

    # design weights
    w = normalize(P[D+1:end])

    # creating a design
    return adjust(Design(x, w))
end


# functions used to convert a design into a particle position
function design2position(ξ::Design)
    # creating a particle position
    P = vcat(ξ.x, ξ.w)

    return P
end


# updating position and velocity of every particle in a swarm
function move!(swarm::Swarm, i::Int64, δ::Float64 = 1e-7)
    # strategic parameters at the current (i-th) iteration
    w  = swarm.W[i]
    c1 = swarm.C1[i]
    c2 = swarm.C2[i]
    
    # current best position in a swarm
    G = swarm.G

    @floop for p in swarm.p
        # a particle's current position
        P = p.P 

        # a particle's best position so far
        B = p.B
        
        r1 = rand(length(P))
        r2 = rand(length(P))
        Vδ = rand(Uniform(-δ, δ), length(P))

        # a particle's current velocity
        V = p.V 

        # calculating a potential
        Φ = [abs(V[d]) + abs(G[d]- P[d]) for d in eachindex(G)]

        # updating particle velocity
        V = [Φ[d] < δ ?  Vδ[d] : w*V[d] + (c1*r1[d])*(B[d] - P[d]) + (c2*r2[d])*(G[d] - P[d]) for d in eachindex(Φ)]

        # updating particle position
        P = P + V
        
        # updating particle information
        setfield!(p, :V, V)
        setfield!(p, :P, P)
    end
end


# adjusting particles' positions after moving according to the constraints 
function adjust_constraints!(swarm::Swarm)
    L = swarm.L
    U = swarm.U

    @floop for p in swarm.p
        P = p.P
        
        # checking constraints and updating position if violated
        D = length(P)
        violated = [k for k in 1:D if (U[k] < P[k]) || (P[k] < L[k])]
        if !isempty(violated)
            for k in violated
                P[k] = U[k] < P[k] ? U[k] : L[k]
            end
            setfield!(p, :P, P)
        end
    end
end


function adjust_mle_constraints!(swarm::Swarm)
    L = swarm.L
    U = swarm.U

    # constraints: Θ={(θ1, θ2, θ3, θ4)∶ θ1 ≥ θ3, θ3<0 and θ2, θ4 >0}
    @floop for p in swarm.p
        θ = p.P
        # checking constraints and updating position if violated
        violated = !all([
            L[2] <= θ[2] <= U[2],
            L[3] <= θ[3] <= U[3],
            L[4] <= θ[4] <= U[4],
            θ[3] <= θ[1] <= U[1]
        ])
        if violated
            # adjusting positions
            θ[2] = L[2] <= θ[2] <= U[2] ? θ[2] : (θ[2] < L[2] ? L[2] : U[2]) # adjusting θ2
            θ[3] = L[3] <= θ[3] <= U[3] ? θ[3] : (θ[3] < L[3] ? L[3] : U[3]) # adjusting θ3
            θ[4] = L[4] <= θ[4] <= U[4] ? θ[4] : (θ[4] < L[4] ? L[4] : U[4]) # adjusting θ4
            θ[1] = θ[3] <= θ[1] <= U[1] ? θ[1] : (θ[1] < θ[3] ? θ[3] : U[1]) # adjusting θ1  
            
            setfield!(p, :P, θ)
        end
    end
end


# evaluating OFV
function evaluate!(ofv::Function, swarm::Swarm)
    @floop for p in swarm.p
        # a particle's current position
        P = p.P

        # evaluating objective function
        f = ofv(P)

        # updating particle's current ofv
        setfield!(p, :f, f)
    end

    # evaluating stopping criteria
    ofvmin, _ = findmin([p.f for p in swarm.p])
    stop = abs(ofvmin - swarm.g) < swarm.ε

    setfield!(swarm, :stop, stop)
end


# updating particles' best position and velocity
function update_best!(swarm::Swarm)
    @floop for p in swarm.p
        # a swarm's best OFV
        g = swarm.g
        
        # a particle's current position
        P = p.P 

        # a particle's current OFV
        f = p.f

        # a particle's best OFV so far
        F = p.F

        # if current OFV is better that the best one so far, update information
        if (f  < F)
            setfield!(p, :F, f)
            setfield!(p, :B, P)
        end

        # if current OFV is better that the global best one so far, update information
        if (f < g)
            setfield!(swarm, :g, f)
            setfield!(swarm, :G, P)
        end
    end
end


# optimizing ofv
function optimize!(ofv::Function, swarm::Swarm)
    N = swarm.N
    i = 1
    while (i <= N && !swarm.stop)
        move!(swarm, i)
        adjust_constraints!(swarm)
        evaluate!(ofv, swarm)
        update_best!(swarm)

        i += 1
    end
    println("number of iterations performed: $(i-1)")

    return position2design(swarm.G)
end


# optimizing likelihood
function optimize_mle!(ofv::Function, swarm::Swarm)
    N = swarm.N
    i = 1
    while (i <= N && !swarm.stop)
        move!(swarm, i)
        adjust_mle_constraints!(swarm)
        evaluate!(ofv, swarm)
        update_best!(swarm)

        i += 1
    end
    println("number of iterations performed: $(i-1)")

    return swarm.G
end


# =============== OFVs for different criteria ===============
# D-optimality    
function ofvD(P::Vector{<:Number}, θ::Vector{<:Number})
    # constructing a design
    ξ = position2design(P)
        
    return ΨD(ξ, θ)
end


# c-optimality for OBD     
function ofvOBD(P::Vector{<:Number}, θ::Vector{<:Number})
    # constructing a design
    ξ = position2design(P)
        
    return Ψobd(ξ, θ)
end
    

# c-optimality for MTD
function ofvMTD(P::Vector{<:Number}, θ::Vector{<:Number}, Γ::Number = 0.2)
    # constructing a design
    ξ = position2design(P)
        
    return Ψmtd(ξ, θ, Γ)
end  