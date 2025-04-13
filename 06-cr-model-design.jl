# ======================================================================
#
# In this script, we define functionality to deal with an experimental 
# designs:
# - a new type `Design` and its constructor(s).
# - `Design` instance representation.
#
# ======================================================================


struct Design
    x::Vector{<:Number}
    w::Vector{<:Number}

    function Design(x::Vector{<:Number}, w::Vector{<:Number})
        @assert length(x) == length(w) "The number of design points must be equat to the number of design weights"
        @assert all(0 .<= w .<= 1) "Weights must be nonnegative numbers."
        @assert isapprox(sum(w), 1.0; atol = 1e-2) "The sum of design weights must be equal to 1."
        
        return new(sort(x), w[sortperm(x)])
    end
end
Design(x::Vector{<:Number}) = Design(x, [1/length(x) for _ in eachindex(x)])
Design(x::Number) = Design([x], [1])


# adjusting a design 
#   - by averaging close points (distance between points ≈ x_distance);
#   - by removing points with low weights (w ≤ w_threshold)
function adjust(ξ::Design; x_distance::Number = 0.1, w_threshold::Number = 0.05, ndigits::Int64 = 5)
    # extracting design points
    x = [round(x, digits = ndigits) for x in ξ.x]

    # number of design points
    nx = length(x)

    # extrcating design weights
    w = ξ.w

    x_adjusted = Number[]
    w_adjusted = Number[]
    
    k = 1
    while (k <= nx)
        xk, wk = x[k], w[k]
        n = 1
        j = k+1
        if (j <= nx)
            while (isapprox(x[k], x[j]; atol = x_distance))
                xk += x[j]
                wk += w[j]
                n += 1
                j += 1
                if (j > nx)
                    break
                end
            end
        end
        push!(w_adjusted, wk)
        push!(x_adjusted, xk/n)
        k = j
    end 

    x_adjusted = [round(x, digits = ndigits) for x in x_adjusted]
    w_adjusted = [round(w, digits = ndigits) for w in w_adjusted]

    x_adjusted = x_adjusted[w_adjusted .> w_threshold]
    w_adjusted = w_adjusted[w_adjusted .> w_threshold]

    return Design(x_adjusted, w_adjusted ./ sum(w_adjusted))
end


# converting an instanse of `Design` into a `String` 
function design2str(ξ::Design)
    x_formatted = [lstrip(@sprintf("% 7.5f", x)) for x in ξ.x]
    w_formatted = [lstrip(@sprintf("% 7.5f", w)) for w in ξ.w]
    
    design = map(arg -> arg[1] * " (" * arg[2] * ")", zip(x_formatted, w_formatted))
    
    return "design (ξ): " * join(design, ", ")
end


# overriding Base function `show` to print out a Design instance
function Base.show(io::IO, ξ::Design)
    println(io, design2str(ξ))
end