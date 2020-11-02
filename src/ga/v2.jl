function ga_v2(
        MODEL_PATH::String,
        objective::Function,
        decode_gene2val::Function,
        nth_param_set::Int64,
        max_generation::Int64,
        n_population::Int64,
        n_children::Int64,
        n_gene::Int64,
        allowable_error::Float64)::Tuple{Array{Float64, 1}, Float64}
    
    #=
        1. Initialization
            As an initial population, create np individuals randomly.
            ga_v2 also represents individuals as n-dimensional real number
            vectors, where n is the dimension of the search space. Set
            Generation to 0, and set the iteration number of converging
            operations Niter to 1.
        2. Selection for reproduction
            As parents for the recombination operator, ENDX, select m
            individuals, p1, p2, . . . ,pm, without replacement from the
            population.
        3. Generation of offsprings
            Generate Nc children by applying ENDX to the selected parents.
            This algorithm assigns the worst objective value to the children.
        4. Local Search (NDM/MGG)
            Apply the local search method to the best individual in a family
            consisting of the two parents, i.e., p1 and p2, and their children.
            Note here that the children are assumed to have the worst objective
            value. Thus, whenever the objective values of the two parents have
            been actually computed in previous generations, the algorithm
            applies the local search to either of the parents. When all of the
            individuals in the family have the same objective value, on the
            other hand, the local search is applied to a randomly selected
            individual from the family.
        5. Selection for survival
            Select two individuals from the family. The first selected
            individual should be the individual with the best objective value,
            and the second should be selected randomly. Then, replace the two
            parents (p1 and p2) with the selected individuals. Note that the
            individual to which the local search has been applied in the
            previous step is always selected as the best.
        6. Application of ENDX/MGG
            To achieve a good search performance, ga_v2 optimizes a function,
            gradually narrowing the search space. For this purpose, the
            converging phase slightly converges the population by repeating the
            following procedure Niter times.
            (i) Select m individuals without replacement from the population.
                The selected individuals, expressed here as p1, p2, . . . , pm,
                are used as the parents for an extended normal distribution
                crossover (ENDX) applied in the next step.
            (ii) Generate Nc children by applying ENDX to the parents selected
                in the previous step. To reduce the computational cost, ga_v2
                forgoes any computation of the objective values of the Nc
                individuals generated here. Instead, the algorithm assigns the
                newly generated children a single objective value, one which is
                inferior to the objective values of any of the possible
                candidate solutions.
            (iii) Select two individuals from a family containing the two
                parents, i.e., p1 and p2, and their children. The first
                selected individual should be the one with the best objective
                value, and the second should be selected randomly. Then,
                replace the two parents with the selected individuals.
        7. Adaptation of Niter
            If the best individual has not improved during the last np
            generations, Niter <- 2 * Niter. Otherwise, set Niter to 1.
        8. Termination
            Stop if the halting criteria are satisfied.
            Otherwise, Generation <- Generation + 1, and return to the step 2.
    =#
    if n_population < n_gene + 2
        error("n_population must be larger than $(n_gene+2)")
    end

    N_iter::Int64 = 1
    N0::Vector{Float64} = zeros(3*n_population)

    population::Matrix{Float64} = get_initial_population(
        MODEL_PATH, objective, nth_param_set, n_population, n_gene
    )
    N0[1] = population[1, end]
    open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "a") do f
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
    while generation <= max_generation
        ip = randperm(n_population)[1:n_gene+2]
        population = converging!(
            objective, ip, population, n_population, n_gene
        )
        population = local_search!(
            objective, ip, population, n_population, n_children, n_gene
        )
        if N_iter > 1
            for _ in 1:N_iter
                ip = randperm(n_population)[1:n_gene+2]
                population = converging!(
                    objective, ip, population, n_population, n_gene
                )
            end
        end

        if generation%length(N0) == 0
            N0[end] = population[1, end]
            if N0[1] == N0[end]
                N_iter *= 2
            else
                N_iter = 1
            end
        else
            N0[generation%length(N0)] = population[1, end]
        end
        open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "a") do f
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


function ga_v2_continue(
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
    if n_population < n_gene + 2
        error("n_population must be larger than $(n_gene+2)")
    end
    
    N_iter::Int64 = 1
    N0::Vector{Float64} = zeros(3*n_population)

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

    N0[1] = population[1, end]
    open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "a") do f
        write(f,
            @sprintf(
                "Generation%d: Best Fitness = %.6e\n",
                count + 1, population[1, end]
            )
        )
    end

    if population[1, end] <= allowable_error
        best_indiv = decode_gene2val(population[1, 1:n_gene])
        best_fitness = population[1, end]

        return best_indiv, best_fitness
    end

    generation::Int64 = 2 + count
    while generation <= max_generation
        ip = randperm(n_population)[1:n_gene+2]
        population = converging!(
            objective, ip, population, n_population, n_gene
        )
        population = local_search!(
            objective, ip, population, n_population, n_children, n_gene
        )
        if N_iter > 1
            for _ in 1:N_iter
                ip = randperm(n_population)[1:n_gene+2]
                population = converging!(
                    objective, ip, population, n_population, n_gene
                )
            end
        end

        if generation%length(N0) == 0
            N0[end] = population[1, end]
            if N0[1] == N0[end]
                N_iter *= 2
            else
                N_iter = 1
            end
        else
            N0[generation%length(N0)] = population[1, end]
        end
        open(strip(MODEL_PATH, '/') * "/logs/$nth_param_set.log", "a") do f
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