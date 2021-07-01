function get_initial_population(
        model::Model,
        nth_param_set::Int64,
        n_population::Int64,
        n_gene::Int64,
        initial_threshold::Float64)::Matrix{Float64}
    open(
        joinpath(
            model.path,
            "logs",
            "$nth_param_set.log"
        ),
        "w",
    ) do f
        write(f, "Generating the initial population...\n\n")
    end
    population::Matrix{Float64} = fill(
        Inf, (n_population, n_gene + 1)
    )
    @inbounds @simd for i = 1:n_population
        while initial_threshold <= population[i, end]
            for j = 1:n_gene
                population[i, j] = rand()
            end
            population[i, end] = model.obj_func(population[i, 1:n_gene])
        end
        open(
            joinpath(
                model.path,
                "logs",
                "$nth_param_set.log"
            ),
            "a",
        ) do f
            write(f, @sprintf("%d / %d\n", i, n_population))
        end
    end
    open(
        joinpath(
            model.path,
            "logs",
            "$nth_param_set.log"
        ),
        "a"
    ) do f
        write(f, "\n----------------------------------------\n")
    end
    population = sortslices(population, dims=1, by=x -> x[end])

    return population
end


function get_initial_population_continue(
        model::Model,
        nth_param_set::Int64,
        n_population::Int64,
        n_gene::Int64,
        initial_threshold::Float64,
        p0_bounds::Vector{Float64})::Matrix{Float64}
    generation::Int64 = readdlm(
        joinpath(
            model.path,
            "fitparam",
            "$nth_param_set",
            "generation.dat"
        )
    )[1, 1]
    best_indiv::Vector{Float64} = readdlm(
        joinpath(
            model.path,
            "fitparam",
            "$nth_param_set",
            "fit_param$generation.dat"
        )
    )[:, 1]
    open(
        joinpath(
            model.path,
            "logs",
            "$nth_param_set.log"
        ),
        "a"
    ) do f
        write(f, 
            "\n----------------------------------------\n" *
            "Generating the initial population...\n"
        )
    end
    population::Matrix{Float64} = fill(
        Inf, (n_population, n_gene + 1)
    )
    @inbounds @simd for i = 1:n_population
        while initial_threshold <= population[i, end]
            for gene_idx = 1:n_gene
                population[i, gene_idx] = model.bestIndivVal2randGene(
                    gene_idx, best_indiv, p0_bounds
                )
                if population[i, gene_idx] > 1.0
                    population[i, gene_idx] = 1.0
                elseif population[i, gene_idx] < 0.0
                    population[i, gene_idx] = 0.0
                end
            end
            population[i, end] = model.obj_func(population[i, 1:n_gene])
        end
        open(
            joinpath(
                model.path,
                "logs",
                "$nth_param_set.log"
            ),
            "a"
        ) do f
            write(f, @sprintf("%d / %d\n", i, n_population))
        end
    end
    open(
        joinpath(
            model.path,
            "logs",
            "$nth_param_set.log"
        ),
        "a"
    ) do f
        write(f, "\n----------------------------------------\n")
    end
    population = sortslices(population, dims=1, by=x -> x[end])

    return population
end