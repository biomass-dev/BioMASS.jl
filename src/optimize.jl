function check_method(local_search_method::String)
    if !(local_search_method in ["mutation", "powell", "de"])
        error(
            "\"$local_search_method\": Invalid local_search_method. "
            * "Should be one of [\"mutaion\", \"Powell\", \"DE\"]"
        )
    elseif local_search_method in ["powell", "de"] && !isinstalled("scipy.optimize")
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
        maxiter::Int64=10,
        local_search_method::String="mutation")
    check_method(
        lowercase(local_search_method)
    )

    for dir in ["fitparam", "logs"]
        if !isdir(
            joinpath(
                model.path,
                dir
            )
        )
            mkdir(
                joinpath(
                    model.path,
                    dir
                )
            )
        end
    end

    try
        files = readdir(
            joinpath(
                model.path,
                "fitparam",
                "$nth_param_set"
            )
        )
        for file in files
            if occursin(".dat", file)
                rm(
                    joinpath(
                        model.path,
                        "fitparam",
                        "$nth_param_set",
                        "$file"
                    )
                )
            end
        end
    catch
        mkdir(
            joinpath(
                model.path,
                "fitparam",
                "$nth_param_set"
            )
        )
    end

    search_rgn::Matrix{Float64} = model.search_region()

    n_population::Int64 = popsize * size(search_rgn, 2)
    n_gene::Int64 = size(search_rgn, 2)

    ga_v2(
        model=model,
        nth_param_set=nth_param_set,
        max_generation=max_generation,
        n_population=n_population,
        n_gene=n_gene,
        allowable_error=allowable_error,
        local_search_method=lowercase(local_search_method),
        n_children=n_children,
        maxiter=maxiter
    )
end


function optimize_continue(
        model::ExecModel,
        nth_param_set::Int64;
        popsize::Int64=5,
        max_generation::Int64=10000,
        allowable_error::Float64=0.0,
        n_children::Int64=50,
        maxiter::Int64=10,
        local_search_method::String="mutation",
        p0_bounds::Vector{Float64}=[0.1, 10.0])
    check_method(
        lowercase(local_search_method)
    )

    search_rgn::Matrix{Float64} = model.search_region()

    n_population::Int64 = popsize * size(search_rgn, 2)
    n_gene::Int64 = size(search_rgn, 2)

    if !isdir(
        joinpath(
            model.path,
            "fitparam",
            "$nth_param_set"
        )
    )
        mkdir(
            joinpath(
                model.path,
                "fitparam",
                "$nth_param_set"
            )
        )

        ga_v2(
            model=model,
            nth_param_set=nth_param_set,
            max_generation=max_generation,
            n_population=n_population,
            n_gene=n_gene,
            allowable_error=allowable_error,
            local_search_method=lowercase(local_search_method),
            n_children=n_children,
            maxiter=maxiter
        )
    else
        ga_v2_continue(
            model=model,
            nth_param_set=nth_param_set,
            max_generation=max_generation,
            n_population=n_population,
            n_gene=n_gene,
            allowable_error=allowable_error,
            local_search_method=lowercase(local_search_method),
            n_children=n_children,
            maxiter=maxiter,
            p0_bounds=p0_bounds
        )
    end
end
