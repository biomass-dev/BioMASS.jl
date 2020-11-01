function ga_v1(
        MODEL_PATH::String,
        objective::Function,
        decode_gene2val::Function,
        nth_param_set::Int64,
        max_generation::Int64,
        n_population::Int64,
        n_children::Int64,
        n_gene::Int64,
        allowable_error::Float64)::Tuple{Array{Float64, 1}, Float64}
    population::Matrix{Float64} = get_initial_population(
        MODEL_PATH, objective, nth_param_set, n_population, n_gene
    )
    open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "w") do f
        write(f,
            @sprintf(
                "Generation%d: Best Fitness = %.6e\n", 1, population[1, end]
            )
        )
    end
    best_indiv::Vector{Float64} = decode_gene2val(population[1, 1:n_gene])
    best_fitness::Float64 = population[1, end]

    open(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/fit_param1.dat", "w") do f
        for val in best_indiv
            write(f, @sprintf("%.6e\n", val))
        end
    end
    open(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/generation.dat", "w") do f
        write(f, @sprintf("%d", 1))
    end
    open(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/best_fitness.dat", "w") do f
        write(f, @sprintf("%.6e", best_fitness))
    end

    if population[1, end] <= allowable_error
        best_indiv = decode_gene2val(population[1, 1:n_gene])
        best_fitness = population[1, end]

        return best_indiv, best_fitness
    end

    generation::Int64 = 2
    while generation < max_generation
        population = mgg_alternation!(
            objective, population, n_population, n_children, n_gene
        )
        open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "w") do f
            write(f,
                @sprintf(
                    "Generation%d: Best Fitness = %.6e\n",
                    generation, population[1, end]
                )
            )
        end
        best_indiv = decode_gene2val(population[1, 1:n_gene])
        if population[1, end] < best_fitness
            open(
                strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/fit_param$generation.dat", "w"
            ) do f
                for val in best_fitness
                    write(f, @sprintf("%.6e\n", val))
                end
            end
            open(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/generation.dat", "w") do f
                write(f, @sprintf("%d", generation))
            end
        end
        best_fitness = population[1, end]
        open(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/best_fitness.dat", "w") do f
            write(f, @sprintf("%.6e", best_fitness))
        end

        if population[1, end] <= allowable_error
            best_indiv = decode_gene2val(population[1, 1:n_gene])
            best_fitness = population[1, end]

            return best_indiv, best_fitness
        end

        open(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/count_num.dat", "w") do f
            write(f, @sprintf("%d", generation))
        end

        generation += 1
    end

    best_indiv = decode_gene2val(population[1, 1:n_gene])
    best_fitness = population[1, end]

    return best_indiv, best_fitness
end


function ga_v1_continue(
        MODEL_PATH::String,
        objective::Function,
        decode_gene2val::Function,
        encode_val2gene::Function,
        encode_bestIndivVal2randGene::Function,
        nth_param_set::Int64,
        max_generation::Int64,
        n_population::Int64,
        n_children::Int64,
        n_gene::Int64,
        allowable_error::Float64, 
        p0_bounds::Vector{Float64})::Tuple{Array{Float64, 1}, Float64}
    count::Int64 = readdlm(
        strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/count_num.dat"
    )[1, 1]
    best_generation::Int64 = readdlm(
        strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/generation.dat"
    )[1, 1]
    best_indiv::Vector{Float64} = readdlm(
        strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/fit_param$best_generation.dat"
    )[:, 1]
    best_indiv_gene::Vector{Float64} = encode_val2gene(best_indiv)
    best_fitness::Float64 = objective(best_indiv_gene)

    population::Matrix{Float64} = get_initial_population_continue(
        MODEL_PATH, objective, encode_bestIndivVal2randGene,
        nth_param_set, n_population, n_gene, p0_bounds
    )
    if best_fitness < population[1, end]
        for i in 1:n_gene
            @inbounds population[1, i] = best_indiv_gene[i]
        end
        population[1, end] = best_fitness
    else
        best_indiv = decode_gene2val(population[1, 1:n_gene])
        best_fitness = population[1, end]
        open(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/fit_param$count.dat", "w") do f
            for i=1:n_gene
                write(f, @sprintf("%.6e", best_indiv[i]))
            end
        end
    end
    open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "w") do f
        write(f,
            @sprintf(
                "Generation%d: Best Fitness = %.6e\n",
                count+1, population[1, end]
            )
        )
    end

    if population[1, end] <= allowable_error
        best_indiv = decode_gene2val(population[1, 1:n_gene])
        best_fitness = population[1, end]

        return best_indiv, best_fitness
    end

    generation::Int64 = 2 + count
    while generation < max_generation
        population = mgg_alternation!(
            objective, population, n_population, n_children, n_gene
        )
        open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "w") do f
            write(f,
                @sprintf(
                    "Generation%d: Best Fitness = %.6e\n",
                    generation, population[1, end]
                )
            )
        end
        best_indiv = decode_gene2val(population[1, 1:n_gene])
        if population[1, end] < best_fitness
            open(
                strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/fit_param$generation.dat"
                , "w"
            ) do f
                for val in best_indiv
                    write(f, @sprintf("%.6e\n", val))
                end
            end
            open(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/generation.dat", "w") do f
                write(f, @sprintf("%d", generation))
            end
        end
        best_fitness = population[1, end]

        open(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/best_fitness.dat", "w") do f
            write(f, @sprintf("%.6e", best_fitness))
        end

        if population[1, end] <= allowable_error
            best_indiv = decode_gene2val(population[1, 1:n_gene])
            best_fitness = population[1, end]

            return best_indiv, best_fitness
        end

        open(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/count_num.dat", "w") do f
            write(f, @sprintf("%d", generation))
        end

        generation += 1
    end

    best_indiv = decode_gene2val(population[1, 1:n_gene])
    best_fitness = population[1, end]

    return best_indiv, best_fitness
end