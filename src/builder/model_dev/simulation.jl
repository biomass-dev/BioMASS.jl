module Sim
include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./set_model.jl")
include("./observable.jl")

using .C
using .V

using Sundials
using SteadyStateDiffEq

# Options for ODE solver
const ABSTOL = 1e-9
const RELTOL = 1e-9
const DTMIN = 1e-8

const normalization = true 
#=
if true, simulation results in each observable 
are divided by their maximum values
=#
const dt = 1.0
t = collect(0.0:dt:100.0)  # 0, 1, 2, ..., 100

const conditions = [

]

simulations = Array{Float64,3}(
    undef, length(observables), length(t), length(conditions)
)


function solveode(
        f::Function,
        u0::Vector{Float64},
        t::Vector{Float64},
        p::Vector{Float64})::Union{ODESolution{}, Nothing}
    local sol::ODESolution{}, is_successful::Bool
    try
        prob = ODEProblem(f,u0,(t[1],t[end]),p)
        sol = solve(
            prob,CVODE_BDF(),
            abstol=ABSTOL,reltol=RELTOL,dtmin=DTMIN,saveat=dt,verbose=false
        )
        is_successful = ifelse(sol.retcode === :Success, true, false)
    catch
        is_successful = false
    finally
        if !is_successful
            GC.gc()
        end
    end
    return is_successful ? sol : nothing
end


function get_steady_state(
        f::Function,
        u0::Vector{Float64},
        p::Vector{Float64})::Vector{Float64}
    local sol::SteadyStateSolution{}, is_successful::Bool
    try
        prob = ODEProblem(diffeq,u0,(0.0,Inf),p)
        prob = SteadyStateProblem(prob)
        sol = solve(
            prob,
            DynamicSS(
                CVODE_BDF();abstol=ABSTOL,reltol=RELTOL
            ),
            dt=dt,dtmin=DTMIN,verbose=false
        )
        is_successful = ifelse(sol.retcode === :Success, true, false)
    catch
        is_successful = false
    finally
        if !is_successful
            GC.gc()
        end
    end
    return is_successful ? sol.u : []
end


function simulate!(p::Vector{Float64}, u0::Vector{Float64})::Union{Bool, Nothing}
    # get steady state
    # 
    # u0 = get_steady_state(diffeq,u0,p)
    # if isempty(u0)
    #    return false
    # end
    for (i,condition) in enumerate(conditions)

        sol = solveode(diffeq,u0,t,p)
        if sol === nothing
            return false
        else
            @inbounds @simd for j in eachindex(t)
                
            end
        end
    end
end
end # module