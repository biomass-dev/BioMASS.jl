function optimize(
        model::ExecModel,
        nth_param_set::Int64;
        n_children::Int64=50,
        max_generation::Int64=10000,
        allowable_error::Float64=0.0)
    
    if !isdir(strip(model.path, '/') * "/fitparam")
        mkdir(strip(model.path, '/') * "/fitparam")
    end
    if !isdir(strip(model.path, '/') * "/logs")
        mkdir(strip(model.path, '/') * "/logs")
    end

    try
        files = readdir(
            strip(model.path, '/') * "/fitparam/$nth_param_set"
        )
        for file in files
            if occursin(".dat",file)
                rm(
                    strip(model.path, '/') * "/fitparam/$nth_param_set/$file"
                )
            end
        end
    catch
        mkdir(strip(model.path, '/') * "/fitparam/$nth_param_set")
    end

    search_rgn::Matrix{Float64} = model.search_region()

    n_population::Int64 = 5*size(search_rgn, 2)
    n_gene::Int64 = size(search_rgn, 2)
    
    ga_v2(
        model,
        nth_param_set,
        max_generation,
        n_population,
        n_children,
        n_gene,
        allowable_error
    )
end


function optimize_continue(
        model::ExecModel,
        nth_param_set::Int64;
        n_children::Int64=50,
        max_generation::Int64=10000,
        allowable_error::Float64=0.0)

    search_rgn::Matrix{Float64} = model.search_region()

    n_population::Int64 = 5*size(search_rgn, 2)
    n_gene::Int64 = size(search_rgn, 2)

    p0_bounds::Vector{Float64} = [0.1, 10.0]  # [lower_bound, upper_bound]

    if !isdir(strip(model.path, '/') * "/fitparam/$nth_param_set")
        mkdir(strip(model.path, '/') * "/fitparam/$nth_param_set")
        
        ga_v2(
            model,
            nth_param_set,
            max_generation,
            n_population,
            n_children,
            n_gene,
            allowable_error
        )
    else
        ga_v2_continue(
            model,
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