function get_initial_population(
        MODEL_PATH::String,
        objective::Function,
        nth_param_set::Int64,
        n_population::Int64,
        n_gene::Int64)::Matrix{Float64}
    open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "w") do f
        write(f, "Generating the initial population...\n\n")
    end
    population::Matrix{Float64} = fill(
        Inf, (n_population, n_gene + 1)
    )
    for i = 1:n_population
        while !isfinite(population[i,end])
            for j = 1:n_gene
                population[i,j] = rand()
            end
            population[i,end] = objective(population[i,1:n_gene])
        end
        open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "a") do f
            write(f, @sprintf("%d / %d\n", i, n_population))
        end
    end
    open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "a") do f
        write(f, "\n----------------------------------------\n")
    end
    population = sortslices(population, dims = 1, by = x->x[end])

    return population
end


function get_initial_population_continue(
        MODEL_PATH::String,
        objective::Function,
        encode_bestIndivVal2randGene::Function,
        nth_param_set::Int64,
        n_population::Int64,
        n_gene::Int64,
        p0_bounds::Vector{Float64})::Matrix{Float64}
    generation::Int64 = readdlm(
        strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/generation.dat"
    )[1,1]
    best_indiv::Vector{Float64} = readdlm(
        strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/fit_param$generation.dat"
    )[:,1]
    open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "a") do f
        write(f,
            "\n----------------------------------------\n"*
            "Generating the initial population...\n"
        )
    end
    population::Matrix{Float64} = fill(
        Inf, (n_population, n_gene + 1)
    )
    for i = 1:n_population
        while !isfinite(population[i,end])
            for j = 1:n_gene
                population[i,j] = encode_bestIndivVal2randGene(
                    j, best_indiv, p0_bounds
                )
                if population[i,j] > 1.0
                    population[i,j] = 1.0
                elseif population[i,j] < 0.0
                    population[i,j] = 0.0
                end
            end
            population[i,end] = objective(population[i,1:n_gene])
        end
        open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "a") do f
            write(f, @sprintf("%d / %d\n", i, n_population))
        end
    end
    open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "a") do f
        write(f, "\n----------------------------------------\n")
    end
    population = sortslices(population, dims = 1, by = x->x[end])

    return population
end