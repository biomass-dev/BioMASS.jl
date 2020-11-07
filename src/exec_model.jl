struct ExecModel
    path::String
    obj_func::Function
    search_region::Function
    gene2val::Function
    val2gene::Function
    bestIndivVal2randGene::Function
end


function set_model_files(dir::Union{String, SubString{String}})
    for child in readdir(dir)
        # if isdir("$dir/$child")
        if child == "name2idx"
            set_model_files("$dir/$child")
        elseif splitext(child)[end] == ".jl"
            include("$dir/$child")
        end
    end
end


function load_model(model_path::String)
    set_model_files(strip(model_path, '/'))
    return ExecModel(
        model_path,
        objective,
        get_search_region,
        decode_gene2val,
        encode_val2gene,
        encode_bestIndivVal2randGene
    )
end