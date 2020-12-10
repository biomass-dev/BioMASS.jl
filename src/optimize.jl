function check_method(local_search_method::String)
    if !(local_search_method in ["mutation", "powell"])
        error(
            "\"$local_search_method\": Invalid local_search_method. "
            * "Should be one of [\"mutaion\", \"powell\"]"
        )
    elseif local_search_method == "powell" && !isinstalled("scipy.optimize")
        error(
            "Cannnot import scipy.optimize functions. Use \"mutation\" for local search."
        )
    end
end


function optimize(
        model::ExecModel,
        nth_param_set::Int64;
        popsize::Int64=5,
        max_generation::Int64=10000,
        allowable_error::Float64=0.0,
        n_children::Int64=50,
        local_search_method::String="mutation")
    check_method(
        lowercase(local_search_method)
    )

    for dir in ["/fitparam", "/logs"]
        if !isdir(strip(model.path, '/') * dir)
            mkdir(strip(model.path, '/') * dir)
        end
    end

    try
        files = readdir(
            strip(model.path, '/') * "/fitparam/$nth_param_set"
        )
        for file in files
            if occursin(".dat", file)
                rm(
                    strip(model.path, '/') * "/fitparam/$nth_param_set/$file"
                )
            end
        end
    catch
        mkdir(strip(model.path, '/') * "/fitparam/$nth_param_set")
    end

    search_rgn::Matrix{Float64} = model.search_region()

    n_population::Int64 = popsize * size(search_rgn, 2)
    n_gene::Int64 = size(search_rgn, 2)

    ga_v2(
        model,
        nth_param_set,
        max_generation,
        n_population,
        n_children,
        n_gene,
        allowable_error,
        lowercase(local_search_method)
    )
end


function optimize_continue(
        model::ExecModel,
        nth_param_set::Int64;
        popsize::Int64=5,
        max_generation::Int64=10000,
        allowable_error::Float64=0.0,
        p0_bounds::Vector{Float64}=[0.1, 10.0],
        n_children::Int64=50,
        local_search_method::String="mutation")
    check_method(
        lowercase(local_search_method)
    )

    search_rgn::Matrix{Float64} = model.search_region()

    n_population::Int64 = popsize * size(search_rgn, 2)
    n_gene::Int64 = size(search_rgn, 2)

    if !isdir(strip(model.path, '/') * "/fitparam/$nth_param_set")
        mkdir(strip(model.path, '/') * "/fitparam/$nth_param_set")

        ga_v2(
            model,
            nth_param_set,
            max_generation,
            n_population,
            n_children,
            n_gene,
            allowable_error,
            lowercase(local_search_method)
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
            p0_bounds,
            lowercase(local_search_method)
        )
    end
end
