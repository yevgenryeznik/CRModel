# ======================================================================
#
# In this script, we define function(s) to visualize CR model
#
# ======================================================================

function plot_cr_model(
    θ::Vector{<:Number}, 
    dmin::Number = log(1e-3),
    dmax::Number = log(1e3);
    Δ::Number = 0.2,
    Γ::Number = 0.2,
    D::Union{Nothing, Vector{<:Number}} = nothing,
    title = "",
    add_bounds = true,
    kwargs...
)
    # calculating MinED, OBD, MTD:
    # on a continuous dose space
    if isnothing(D) 
        MinED = round(mined(θ, Δ)[1], digits = 2)
        OBD = round(obd(θ)[1], digits = 2)
        MTD = round(mtd(θ, Γ)[1], digits = 2)
        D1 = [MinED, OBD, MTD]
    # on a discrete dose space      
    else 
        MinED = mined(θ, D, Δ)
        OBD = obd(θ, D)
        MTD = mtd(θ, D, Γ)
        D1 = D                  
    end

    d = collect(LinRange(dmin, dmax, 50))
    yd = hcat(η(d, θ)...)
    yD1 = hcat(η(D1, θ)...)
    


    lbl = [latexify("π_$k(d, θ)") for k in [0, 1, 2]]
    clr = [:blue, :green, :red]
    mkr = [:circle, :rect, :utriangle]

    xticks_value = [round(dmin, digits = 3); D1; round(dmax, digits = 3)]
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

    θstr = map(i -> "θ_$i=$(θ[i])", eachindex(θ))

    fig = Figure(
        fontsize = 18,
        size = (1200, 800),
        dpi = 300
    )
    ax = Axis(
        fig[1, 1], 
        title = latexstring("\\textbf{$title} (" * join(θstr, ",\\:") * ")"),
        titlesize = 24,
        xlabel = L"$\log{\left(\text{dose}\right)}$",
        ylabel = "probability",
        xticks = (xticks_value, xticks_label),
        yticks = 0:0.1:1,
        xticksize = 0,
        yticksize = 0,
        xticklabelrotation = π/6,
        kwargs...
    )

    for r in axes(yd, 1)
        lines!(ax, d, yd[r, :], color = clr[r], label = lbl[r])
        scatter!(ax, D1, yD1[r, :], color = clr[r], marker = mkr[r], label = lbl[r], markersize = 20)
    end
    axislegend(ax, "", merge = true, framevisible = false, labelsize = 20)
    hidespines!(ax)

    if add_bounds
        vlines!(ax, [D1[1], D1[end]], color = :black)
    end

    return fig, ax    
end