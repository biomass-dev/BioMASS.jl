module ExecModel

using Printf
using LinearAlgebra

export
    C,
    V,
    observables,
    observables_index,
    Sim,
    Exp,
    param_values,
    initial_values,
    objective,
    get_search_index,
    get_search_region,
    decode_gene2val,
    encode_val2gene,
    encode_bestIndivVal2randGene,
    update_param

include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./set_model.jl")
include("./observable.jl")
include("./experimental_data.jl")
include("./simulation.jl")
include("./fitness.jl")
include("./set_search_param.jl")

end  # module