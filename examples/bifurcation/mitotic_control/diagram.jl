include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./set_model.jl")
include("./forwarddiff.jl")

const BP = C.InhCDK     # name(index) of bifurcation parameter

const SN = V.NUM        # num of state variables
const PN = 1            # num of parameters
const VN = SN + PN      # num of variables

function calc_fixed_point_vec(model_path::String)
    global fp, br

    fp = []
    ev = []
    br = []

    for i in 1:4
        global p = param_values()
        # i==1 -> Control
        if i == 2 # Wee1 inhibition
            p[C.kweeS] = 0.0
            p[C.kweeF] = 0.0
        elseif i == 3 # Gwl siRNA
            p[C.Gwtot] = 0.0
        elseif i == 4 # Wee1 inhibition & Gwl siRNA
            p[C.kweeS] = 0.0
            p[C.kweeF] = 0.0
            p[C.Gwtot] = 0.0
        end

        new_curve!(
            model_path, p, diffeq2, get_derivatives, get_steady_state,
            direction=false, bifparam=BP, n_state=SN
        )
        push!(fp,readdlm(model_path * "/data/fp.dat",'\t',Float64,'\n'))
        push!(ev,readdlm(model_path * "/data/ev.dat",'\t',Float64,'\n'))
        push!(br,get_bistable_regime(ev[i], SN))
    end
end


function bifurcation_diagram(model_path::String)
    rc("figure",figsize = (20,3))
    rc("font",family = "Arial")
    rc("font",size = 14)
    rc("axes",linewidth = 1.2)
    rc("xtick.major",width = 1.2)
    rc("ytick.major",width = 1.2)
    rc("lines",linewidth = 2)

    for (i,(fixed_point,unstable_ss)) in enumerate(zip(fp,br))
        if length(unstable_ss) > 0
            intermediate_ss = []
            for j=2:length(unstable_ss)
                if unstable_ss[j] - unstable_ss[j-1] != 1
                    intermediate_ss = unstable_ss[j-1]:unstable_ss[j]
                end
            end
        end

        subplot(1,4,i)
        if i == 1
            plot(
                fixed_point[1:unstable_ss[1]-1,VN+1],
                fixed_point[1:unstable_ss[1]-1,V.Subp+1],
                color="royalblue"
            )
            plot(
                fixed_point[unstable_ss[1]:intermediate_ss[1]-1,VN+1],
                fixed_point[unstable_ss[1]:intermediate_ss[1]-1,V.Subp+1],
                color="darkgray","--"
            )
            plot(
                fixed_point[intermediate_ss,VN+1],
                fixed_point[intermediate_ss,V.Subp+1],
                color="darkorange"
            )
            plot(
                fixed_point[intermediate_ss[end]+1:unstable_ss[end],VN+1],
                fixed_point[intermediate_ss[end]+1:unstable_ss[end],V.Subp+1],
                color="darkgray","--"
            )
            plot(
                fixed_point[unstable_ss[end]+1:end,VN+1],
                fixed_point[unstable_ss[end]+1:end,V.Subp+1],
                color="crimson"
            )
            xlim(0,0.81)
            xlabel("1NMPP1 (μM)")
            ylim(-0.05,1.05)
            yticks([0,0.5,1],[0,50,100])
            ylabel("Sub-p (%)")
            title("Control",fontsize=18)
        elseif i == 2
            plot(
                fixed_point[1:unstable_ss[1]-1,VN+1],
                fixed_point[1:unstable_ss[1]-1,V.Subp+1],
                color="royalblue"
            )
            plot(
                fixed_point[unstable_ss,VN+1],
                fixed_point[unstable_ss,V.Subp+1],
                color="darkgray","--"
            )
            plot(
                fixed_point[unstable_ss[end]+1:end,VN+1],
                fixed_point[unstable_ss[end]+1:end,V.Subp+1],
                color="crimson"
            )
            xlim(0,0.81)
            xlabel("1NMPP1 (μM)")
            ylim(-0.05,1.05)
            yticks([0,0.5,1],[0,50,100])
            title("Wee1 inhibition",fontsize=18)
        elseif i == 3
            plot(
                fixed_point[1:unstable_ss[1]-1,VN+1],
                fixed_point[1:unstable_ss[1]-1,V.Subp+1],
                color="royalblue"
            )
            plot(
                fixed_point[unstable_ss,VN+1],
                fixed_point[unstable_ss,V.Subp+1],
                color="darkgray","--"
            )
            plot(
                fixed_point[unstable_ss[end]+1:end,VN+1],
                fixed_point[unstable_ss[end]+1:end,V.Subp+1],
                color="crimson"
            )
            xlim(0,0.81)
            xlabel("1NMPP1 (μM)")
            ylim(-0.05,1.05)
            yticks([0,0.5,1],[0,50,100])
            title("Gwl siRNA",fontsize=18)
        elseif i == 4
            plot(
                fixed_point[:,VN+1],
                fixed_point[:,V.Subp+1],
                color="darkorange"
            )
            xlim(0,0.81)
            xlabel("1NMPP1 (μM)")
            ylim(-0.05,1.05)
            yticks([0,0.5,1],[0,50,100])
            title("Wee1 inhibition & Gwl siRNA",fontsize=18)
        end
    end
    savefig(model_path * "/bifurcation_diagram.pdf", bbox_inches="tight")
    close()
end