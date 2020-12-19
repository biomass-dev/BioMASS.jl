#= 
Yao, G., Lee, T. J., Mori, S., Nevins, J. R. & You, L. A bistable Rb-E2F switch 
underlies the restriction point. Nat. Cell Biol. 10, 476–482 (2008). 
https://doi.org/10.1038/ncb1711 =#

include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./set_model.jl")
include("./forwarddiff.jl")

const BP = C.S      # name(index) of bifurcation parameter

const SN = V.NUM    # num of state variables
const PN = 1        # num of parameters
const VN = SN + PN  # num of variables


function calc_fixed_point_vec(model_path::String)::Tuple{Array,Array}
    global p = param_values()
    new_curve!(
        model_path, p, diffeq2, get_derivatives, get_steady_state,
        direction=false, bifparam=BP, n_state=SN
    )
    fp::Array = readdlm(model_path * "/data/fp.dat", '\t', Float64, '\n')
    ev::Array = readdlm(model_path * "/data/ev.dat", '\t', Float64, '\n')
    br::Array = get_bistable_regime(ev, SN)

    return fp, br
end

function bifurcation_diagram(model_path::String, fp::Array, br::Array)
    rc("figure", figsize=(8, 6))
    rc("font", family="Arial")
    rc("font", size=24)
    rc("axes", linewidth=1)
    rc("xtick.major", width=1)
    rc("ytick.major", width=1)
    rc("lines", linewidth=3)

    plot(fp[1:br[1] - 1,VN + 1], fp[1:br[1] - 1,V.E + 1], "k-")
    plot(fp[br,VN + 1], fp[br,V.E + 1], lw=1.5, "k--")
    plot(fp[br[end] + 1:end,VN + 1], fp[br[end] + 1:end,V.E + 1], "k-")

    xlabel("Serum (percentage)")
    ylabel("E2F (μM)")

    xlim(0, 2)
    xticks([0,0.5,1,1.5,2])
    yscale("log")
    ylim(1e-4, 2)
    yticks([1e-4,1e-2,1])

    savefig(model_path * "/bifurcation_diagram.pdf", bbox_inches="tight")
    close()
end