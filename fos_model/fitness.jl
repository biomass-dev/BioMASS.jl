# Residual Sum of Squares
function compute_objval_rss(
        sim_data::Vector{Float64},
        exp_data::Vector{Float64})::Float64
    error::Float64 = 0.0
    for i in eachindex(exp_data)
        @inbounds error += (sim_data[i] - exp_data[i])^2
    end
    return error
end


# Cosine similarity
function compute_objval_cos(
        sim_data::Vector{Float64},
        exp_data::Vector{Float64})::Float64
    error::Float64 = 1.0 - dot(sim_data,exp_data)/(norm(sim_data)*norm(exp_data))
    return error
end


function conditions_index(condition_name::String)::Int
    return findfirst(isequal(condition_name), Sim.conditions)
end


function diff_sim_and_exp(
        sim_matrix::Matrix{Float64},
        exp_dict::Dict{String,Array{Float64,1}},
        exp_timepoint::Vector{Float64},
        conditions::Vector{String};
        sim_norm_max::Float64)::Tuple{Vector{Float64}, Vector{Float64}}
    sim_result::Vector{Float64} = []
    exp_result::Vector{Float64} = []

    for (idx,condition) in enumerate(conditions)
        if condition in keys(exp_dict)
            append!(sim_result,sim_matrix[Int.(exp_timepoint.+1),idx])
            append!(exp_result,exp_dict[condition])
        end
    end

    return (sim_result./sim_norm_max, exp_result)
end


# Define an objective function to be minimized.
function objective(indiv_gene::Vector{Float64})::Float64    
    indiv::Vector{Float64} = decode_gene2val(indiv_gene)

    (p,u0) = update_param(indiv)

    if Sim.simulate!(p,u0) isa Nothing
        error::Vector{Float64} = zeros(length(observables))
        for (i,obs_name) in enumerate(observables)
            if isassigned(Exp.experiments,i)
                if length(Sim.normalization) > 0
                    norm_max::Float64 = (
                        Sim.normalization[obs_name]["timepoint"] !== nothing ? maximum(
                            Sim.simulations[
                                i,
                                Sim.normalization[obs_name]["timepoint"],
                                [conditions_index(c) for c in Sim.normalization[obs_name]["condition"]]
                            ]
                        ) : maximum(
                            Sim.simulations[
                                i,
                                :,
                                [conditions_index(c) for c in Sim.normalization[obs_name]["condition"]]
                            ]
                        )
                    )
                end
                error[i] = compute_objval_rss(
                    diff_sim_and_exp(
                        Sim.simulations[i,:,:],
                        Exp.experiments[i],
                        Exp.get_timepoint(obs_name),
                        Sim.conditions,
                        sim_norm_max = (
                            length(Sim.normalization) == 0 ? 1.0 : norm_max
                        )
                    )...
                )
            end
        end
        #=
        error::Vector{Float64} = zeros(16)

        sim_norm_max = maximum(Sim.simulations[observables_index("Phosphorylated_MEKc"),:,:])
        exp_t = Exp.get_timepoint(observables_index("Phosphorylated_MEKc"))
        error[1] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_MEKc"), Int.(exp_t.+1), conditions_index("EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_MEKc")]["EGF"]
        )
        error[2] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_MEKc"), Int.(exp_t.+1), conditions_index("HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_MEKc")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[observables_index("Phosphorylated_ERKc"),:,:])
        exp_t = Exp.get_timepoint(observables_index("Phosphorylated_ERKc"))
        error[3] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_ERKc"), Int.(exp_t.+1), conditions_index("EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_ERKc")]["EGF"]
        )
        error[4] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_ERKc"), Int.(exp_t.+1), conditions_index("HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_ERKc")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[observables_index("Phosphorylated_RSKw"),:,:])
        exp_t = Exp.get_timepoint(observables_index("Phosphorylated_RSKw"))
        error[5] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_RSKw"), Int.(exp_t.+1), conditions_index("EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_RSKw")]["EGF"]
        )
        error[6] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_RSKw"), Int.(exp_t.+1), conditions_index("HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_RSKw")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[observables_index("Phosphorylated_CREBw"),:,:])
        exp_t = Exp.get_timepoint(observables_index("Phosphorylated_CREBw"))
        error[7] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_CREBw"), Int.(exp_t.+1), conditions_index("EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_CREBw")]["EGF"]
        )
        error[8] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_CREBw"), Int.(exp_t.+1), conditions_index("HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_CREBw")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[observables_index("dusp_mRNA"),:,:])
        exp_t = Exp.get_timepoint(observables_index("dusp_mRNA"))
        error[9] = compute_objval_rss(
            Sim.simulations[observables_index("dusp_mRNA"), Int.(exp_t.+1), conditions_index("EGF")]./sim_norm_max,
            Exp.experiments[observables_index("dusp_mRNA")]["EGF"]
        )
        error[10] = compute_objval_rss(
            Sim.simulations[observables_index("dusp_mRNA"), Int.(exp_t.+1), conditions_index("HRG")]./sim_norm_max,
            Exp.experiments[observables_index("dusp_mRNA")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[observables_index("cfos_mRNA"),:,:])
        exp_t = Exp.get_timepoint(observables_index("cfos_mRNA"))
        error[11] = compute_objval_rss(
            Sim.simulations[observables_index("cfos_mRNA"), Int.(exp_t.+1), conditions_index("EGF")]./sim_norm_max,
            Exp.experiments[observables_index("cfos_mRNA")]["EGF"]
        )
        error[12] = compute_objval_rss(
            Sim.simulations[observables_index("cfos_mRNA"), Int.(exp_t.+1), conditions_index("HRG")]./sim_norm_max,
            Exp.experiments[observables_index("cfos_mRNA")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[observables_index("cFos_Protein"),:,:])
        exp_t = Exp.get_timepoint(observables_index("cFos_Protein"))
        error[13] = compute_objval_rss(
            Sim.simulations[observables_index("cFos_Protein"), Int.(exp_t.+1), conditions_index("EGF")]./sim_norm_max,
            Exp.experiments[observables_index("cFos_Protein")]["EGF"]
        )
        error[14] = compute_objval_rss(
            Sim.simulations[observables_index("cFos_Protein"), Int.(exp_t.+1), conditions_index("HRG")]./sim_norm_max,
            Exp.experiments[observables_index("cFos_Protein")]["HRG"]
        )

        sim_norm_max = maximum(Sim.simulations[observables_index("Phosphorylated_cFos"),:,:])
        exp_t = Exp.get_timepoint(observables_index("Phosphorylated_cFos"))
        error[15] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_cFos"), Int.(exp_t.+1), conditions_index("EGF")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_cFos")]["EGF"]
        )
        error[16] = compute_objval_rss(
            Sim.simulations[observables_index("Phosphorylated_cFos"), Int.(exp_t.+1), conditions_index("HRG")]./sim_norm_max,
            Exp.experiments[observables_index("Phosphorylated_cFos")]["HRG"]
        )
        =#
        return sum(error)
    else
        return Inf
    end
end