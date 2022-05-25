module Sim
include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./set_model.jl")
include("./observable.jl")

using .C
using .V

using DelayDiffEq

normalization = Dict{String,Dict{}}()
for observable in observables
    normalization[observable] = Dict(
        "timepoint" => nothing,
        "condition" => ["WT"]
    )
end

const dt = 1.0

t = collect(0.0:1.0:360.0)  # 0, 1, 2, ..., 360 [min.]

const sstime = 1000.0  # time to reach steady state

const conditions = ["WT"]

simulations = Array{Float64,3}(
    undef, length(observables), length(conditions), length(t)
)

function solvedde(
    f::Function, u0::Vector{Float64}, history::Vector{Float64},
    tspan::Tuple{Float64,Float64}, p::Vector{Float64}, tau::Float64)
    h(p, t) = history
    lags = [tau]
    prob = DDEProblem(f, u0, h, tspan, p; constant_lags=lags)
    alg = MethodOfSteps(BS3())
    sol = solve(
        prob, alg, saveat=dt, abstol=1e-8, reltol=1e-8, verbose=false
    )
    return sol
end


function get_steady_state(
    p::Vector{Float64}, u0::Vector{Float64},
    sstime::Float64, tau::Float64)::Vector{Float64}
    # get steady state (t<0)
    p[C.term] = 1.0
    history::Vector{Float64} = u0
    tspan::Tuple{Float64,Float64} = (0.0, sstime)
    try
        sol = solvedde(diffeq!, u0, history, tspan, p, tau)
        if sol.retcode === :Success
            return sol[:, end]
        else
            return []
        end
    catch
        return []
    end
end


function get_time_course(
    p::Vector{Float64}, u0::Vector{Float64},
    sstime::Float64, tau::Float64)
    p1::Vector{Float64} = copy(p)
    p1[C.term] = 0.0
    u1::Vector{Float64} = get_steady_state(p, u0, sstime, tau)
    if isempty(u1)
        return nothing
    end
    history::Vector{Float64} = u1
    tspan::Tuple{Float64,Float64} = (0.0, t[end])
    try
        sol = solvedde(diffeq!, u1, history, tspan, p1, tau)
        return ifelse(sol.retcode === :Success, sol, nothing)
    catch
        return nothing
    end
end


function simulate!(
    p::Vector{Float64},
    u0::Vector{Float64})::Union{Bool,Nothing}
    for (i, condition) in enumerate(conditions)
        # if condition == "WT"
        #    pass
        # end
        sol = get_time_course(p, u0, sstime, p[C.delayrnae])
        if sol === nothing
            return false
        else
            @inbounds @simd for j in eachindex(t)
                simulations[observables_index("Nuclear_NFkB"), i, j] = (
                    sol.u[i][V.NFKBn]
                )
            end
        end
    end
end
end # module