struct ExecModel
    path::String
    parameters::Module
    species::Module
    obj_func::Function
    search_idx::Function
    search_region::Function
    gene2val::Function
    val2gene::Function
    bestIndivVal2randGene::Function
end
function set_model_files(dir::String)
    for child in readdir(dir)
        # if isdir("$dir/$child")
        if child == "name2idx"
            set_model_files(
                joinpath(
                    "$dir", "$child"
                )
            )
        elseif splitext(child)[end] == ".jl"
            include(
                joinpath(
                    "$dir", "$child"
                )
            )
        end
    end
end
function ExecModel(model_path::String)
    set_model_files(model_path)
    ExecModel(
        model_path,
        C,
        V,
        objective,
        get_search_index,
        get_search_region,
        decode_gene2val,
        encode_val2gene,
        encode_bestIndivVal2randGene
    )
end

function load_model(path_to_model::String)
    if isdir(path_to_model)
        return ExecModel(path_to_model)
    else
        error("$path_to_model: No such directory")
    end
end