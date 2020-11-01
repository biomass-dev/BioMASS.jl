module Exp
include("./observable.jl")

experiments = Array{Dict{String,Array{Float64,1}},1}(undef, length(observables))
error_bars = Array{Dict{String,Array{Float64,1}},1}(undef, length(observables))


function get_timepoint(obs_name::String)::Vector{Float64}
    return []
end
end # module