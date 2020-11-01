import .Main

function optimize(
        MODEL_PATH::String,
        nth_param_set::Int64;
        max_generation::Int64=10000,
        allowable_error::Float64=0.0)
    
    if !isdir(strip(MODEL_PATH, '/') * "/fitparam")
        mkdir(strip(MODEL_PATH, '/') * "/fitparam")
    end
    if !isdir(strip(MODEL_PATH, '/') * "/logs")
        mkdir(strip(MODEL_PATH, '/') * "/logs")
    end

    try
        files = readdir(
            strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set"
        )
        for file in files
            if occursin(".dat",file)
                rm(
                    strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set/$file"
                )
            end
        end
    catch
        mkdir(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set")
    end

    search_rgn::Matrix{Float64} = Main.ExecModel.get_search_region()

    n_population::Int64 = 5*size(search_rgn, 2)
    n_children::Int64 = 50
    n_gene::Int64 = size(search_rgn, 2)
    
    (best_indiv, best_fitness) = ga_v2(
        MODEL_PATH,
        Main.ExecModel.objective,
        Main.decode_gene2val,
        nth_param_set,
        max_generation,
        n_population,
        n_children,
        n_gene,
        allowable_error
    )
end


function optimize_continue(
        MODEL_PATH::String,
        nth_param_set::Int64;
        max_generation::Int64=10000,
        allowable_error::Float64=0.0)

    search_rgn::Matrix{Float64} = Main.ExecModel.get_search_region()

    n_population::Int64 = 5*size(search_rgn, 2)
    n_children::Int64 = 50
    n_gene::Int64 = size(search_rgn, 2)

    p0_bounds::Vector{Float64} = [0.1, 10.0]  # [lower_bound, upper_bound]

    if !isdir(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set")
        mkdir(strip(MODEL_PATH, '/') * "/fitparam/$nth_param_set")
        
        (best_indiv, best_fitness) = ga_v2(
            MODEL_PATH,
            Main.ExecModel.objective,
            Main.decode_gene2val,
            nth_param_set,
            max_generation,
            n_population,
            n_children,
            n_gene,
            allowable_error
        )
    else
        (best_indiv, best_fitness) = ga_v2_continue(
            MODEL_PATH,
            Main.ExecModel.objective,
            Main.decode_gene2val,
            Main.encode_val2gene,
            Main.encode_bestIndivVal2randGene,
            nth_param_set,
            max_generation,
            n_population,
            n_children,
            n_gene,
            allowable_error,
            p0_bounds
        )
    end
end