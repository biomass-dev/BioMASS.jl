module BioMASS

using Printf
using LinearAlgebra
using StatsBase
using Statistics
using DelimitedFiles

export
    # Parameter Estimation
    optimize,
    optimize_continue,
    param2biomass,
    ExecModel,
    load_model,
    visualize,
    # Bifurcation Analysis
    create_diffeq,
    new_curve!,
    get_bistable_regime

function isinstalled(pymodule::String)::Bool
    try
        pyimport(pymodule)
        return true
    catch
        return false
    end
end

include("exec_model.jl")
include("convert.jl")
include("optimize.jl")
include("ga/initial_population.jl")
include("ga/converging.jl")
include("ga/local_search.jl")
include("ga/v2.jl")
if isinstalled("matplotlib")
    include("visualize.jl")
else
    function visualize(model::ExecModel; kwargs...)
        error(
            "The Python package matplotlib could not be imported by pyimport.\n"
            * "Usually this means that you did not install matplotlib in the "
            * "Python version being used by PyCall."
        )
    end
end
include("continuation.jl")
end # module
