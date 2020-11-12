module BioMASS

using Printf
using LinearAlgebra
using Random
using StatsBase
using Statistics
using DelimitedFiles

export
    optimize,
    optimize_continue,
    param2biomass,
    ExecModel,
    load_model,
    isinstalled_plt

function isinstalled_plt()::Bool
    try
        pyimport("matplotlib")
        return true
    catch
        return false
    end
end

include("exec_model.jl")
include("convert.jl")
include("optimize.jl")
include("ga/initial_population.jl")
include("ga/undxmgg.jl")
include("ga/converging.jl")
include("ga/local_search.jl")
include("ga/v1.jl")
include("ga/v2.jl")
if isinstalled_plt()
    include("visualize.jl")
    export visualize
end
end
