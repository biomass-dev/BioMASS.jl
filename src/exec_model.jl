struct ExecModel
    path::String
    obj_func::Function
    search_region::Function
    gene2val::Function
    val2gene::Function
    bestIndivVal2randGene::Function
end

function load_model(model_path::String)
    include(strip(model_path, '/') * "/name2idx/parameters.jl")
    include(strip(model_path, '/') * "/name2idx/species.jl")
    include(strip(model_path, '/') * "/set_model.jl")
    include(strip(model_path, '/') * "/observable.jl")
    include(strip(model_path, '/') * "/experimental_data.jl")
    include(strip(model_path, '/') * "/simulation.jl")
    include(strip(model_path, '/') * "/fitness.jl")
    include(strip(model_path, '/') * "/set_search_param.jl")

    return ExecModel(
        model_path,
        objective,
        get_search_region,
        decode_gene2val,
        encode_val2gene,
        encode_bestIndivVal2randGene
    )
end