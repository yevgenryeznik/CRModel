function check_optimality(
    ξ::Design, 
    θ::Vector{<:Number}, 
    ψ::Function, 
    X::Vector{<:Number},
    MinED::Number,
    OBD::Number,
    MTD::Number;
    title::String = "",
    show_title::Bool = false,
    show_xlabel::Bool = false,
    show_legend::Bool = false
)
    θstr = map(i -> latexstring("\\theta_$i=$(θ[i])"), eachindex(θ))
    lb, ub = first(X), last(X)    
    xval = collect(LinRange(lb, ub, 100))
    yval = ψ.(xval)

    digs = 3
    x = round.(ξ.x, digits = digs)
    xticks_value = unique(sort([X; MinED; OBD; MTD; x]))
    xticks_label = ["$v" for v in xticks_value]

    for i in eachindex(xticks_value)
        if xticks_value[i] == MinED
            xticks_label[i] = xticks_label[i]  * "\n(MinED)"
        end
        if xticks_value[i] == OBD
            xticks_label[i] = xticks_label[i]  * "\n(OBD)"
        end
        if xticks_value[i] == MTD
            xticks_label[i] = xticks_label[i]  * "\n(MTD)"
        end
    end

    fig = Figure(
        fontsize = 18,
        size = (1200, 800),
        dpi = 300   
    )
    ax = Axis(
        fig[1, 1], 
        title = show_title ? latexstring("\\textbf{$title} (" * join(θstr, ",\\:") * ")") : "",
        titlesize = 24,
        xlabel = show_xlabel ? L"$d = \log{\left(\text{dose}\right)}$" : "",
        ylabel = "",
        xticks = (xticks_value, xticks_label),
        xticksize = 0,
        yticksize = 0,
        xticklabelrotation = π/6,
        limits = (nothing, nothing, nothing, 0.1)
    )
    lines!(ax, xval, yval, label = "sensitivity function", linewidth = 5)
    scatter!(ax, x, zeros(length(x)), color = :green, markersize = 40, label = "design points")
    hlines!(ax, [0], linestyle = :dash, linewidth = 5, color = :black)
    hidespines!(ax)
    if show_legend
        axislegend(ax, "", merge = true, framevisible = false, labelsize = 30, position = :rb)
    end

    # adding weights inforamtion
    w = round.(ξ.w, digits = digs)
    for (xi, wi) in zip(x, w)
        text!(ax, xi, 0.01, text = string(wi), fontsize = 30, rotation = π/6)
    end
    
    fig, ax
end