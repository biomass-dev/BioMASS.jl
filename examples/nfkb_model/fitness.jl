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
    error::Float64 = 1.0 - dot(sim_data, exp_data) / (norm(sim_data) * norm(exp_data))
    return error
end


function conditions_index(condition_name::String)::Int
    if !(condition_name in Sim.conditions)
        error("$condition_name is not defined in Sim.conditions")
    end
    return findfirst(isequal(condition_name), Sim.conditions)
end


function diff_sim_and_exp(
    sim_matrix::Matrix{Float64},
    exp_dict::Dict{String,Array{Float64,1}},
    exp_timepoint::Vector{Float64},
    conditions::Vector{String};
    sim_norm_max::Float64)::Tuple{Vector{Float64},Vector{Float64}}
    sim_result::Vector{Float64} = []
    exp_result::Vector{Float64} = []

    for (idx, condition) in enumerate(conditions)
        if condition in keys(exp_dict)
            append!(sim_result, sim_matrix[Int.(exp_timepoint .+ 1),idx])
            append!(exp_result, exp_dict[condition])
        end
    end

    return (sim_result ./ sim_norm_max, exp_result)
end


# Define an objective function to be minimized.
function objective(indiv_gene)::Float64    
    indiv::Vector{Float64} = decode_gene2val(indiv_gene)

    (p, u0) = update_param(indiv)

    if Sim.simulate!(p, u0) isa Nothing
        error::Vector{Float64} = zeros(length(observables))
        for (i, obs_name) in enumerate(observables)
            if isassigned(Exp.experiments, i)
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
                    sim_norm_max=ifelse(
                        length(Sim.normalization) == 0, 1.0, norm_max
                    )
                )...
            )
            end
        end
        return sum(error) # < 1e12
    else
        return 1e12
    end
end