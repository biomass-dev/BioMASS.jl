module Sim
include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./ode.jl")
include("./observable.jl")

using .C
using .V

using Sundials
using SteadyStateDiffEq

# Options for ODE solver
const ABSTOL = 1e-8
const RELTOL = 1e-8

normalization = Dict{String,Dict{}}()
for observable in observables
    normalization[observable] = Dict(
        "timepoint" => nothing,
        "condition" => ["EGF", "HRG"]
    )
end

const dt = 1.0
const t = collect(0.0:dt:5400.0)  # 0, 1, 2, ..., 5400 [sec.]

const conditions = ["EGF", "HRG"]

simulations = Array{Float64,3}(
    undef, length(observables), length(conditions), length(t)
)


function solveode(
    f::Function,
    u0::Vector{Float64},
    t::Vector{Float64},
    p::Vector{Float64})::Union{ODESolution{},Nothing}
    local sol::ODESolution{}, is_successful::Bool
    prob = ODEProblem(f, u0, (t[1], t[end]), p)
    try
        prob = ODEProblem(f, u0, (t[1], t[end]), p)
        sol = solve(
            prob, CVODE_BDF(),
            abstol=ABSTOL,
            reltol=RELTOL,
            saveat=dt,
            dtmin=eps(),
            verbose=false
        )
        #is_successful = ifelse(sol.retcode === :Success, true, false)
        return sol
    catch
        #is_successful = false
        GC.gc()
        #finally
        #    if !is_successful
        #        GC.gc()
        #    end
        return nothing
    end
    #return is_successful ? sol : nothing
end


function get_steady_state(
    f::Function,
    u0::Vector{Float64},
    p::Vector{Float64})::Vector{Float64}
    local sol::SteadyStateSolution{}, is_successful::Bool
    try
        prob = ODEProblem(f, u0, (0.0, Inf), p)
        prob = SteadyStateProblem(prob)
        sol = solve(
            prob,
            DynamicSS(
                CVODE_BDF();
                abstol=ABSTOL,
                reltol=RELTOL
            ),
            dt=dt,
            dtmin=eps(),
            verbose=false
        )
        #is_successful = ifelse(sol.retcode === :Success, true, false)
        return sol.u
    catch
        #is_successful = false
        GC.gc()
        #finally
        #    if !is_successful
        #        GC.gc()
        #    end
        return []
    end
    #return is_successful ? sol.u : []
end


function simulate!(p::Vector{Float64}, u0::Vector{Float64})::Union{Bool,Nothing}
    # get steady state
    p[C.Ligand] = p[C.no_ligand]
    u0 = get_steady_state(diffeq!, u0, p)
    if isempty(u0)
        return false
    end
    # add ligand
    for (i, condition) in enumerate(conditions)
        if condition == "EGF"
            p[C.Ligand] = p[C.EGF]
        elseif condition == "HRG"
            p[C.Ligand] = p[C.HRG]
        end
        sol = solveode(diffeq!, u0, t, p)
        if sol === nothing
            return false
        else
            @inbounds @simd for j in eachindex(t)
                simulations[observables_index("Phosphorylated_MEKc"), i, j] = (
                    sol.u[j][V.ppMEKc]
                )
                simulations[observables_index("Phosphorylated_ERKc"), i, j] = (
                    sol.u[j][V.pERKc] + sol.u[j][V.ppERKc]
                )
                simulations[observables_index("Phosphorylated_RSKw"), i, j] = (
                    sol.u[j][V.pRSKc] + sol.u[j][V.pRSKn] * (p[C.Vn] / p[C.Vc])
                )
                simulations[observables_index("Phosphorylated_CREBw"), i, j] = (
                    sol.u[j][V.pCREBn] * (p[C.Vn] / p[C.Vc])
                )
                simulations[observables_index("dusp_mRNA"), i, j] = (
                    sol.u[j][V.duspmRNAc]
                )
                simulations[observables_index("cfos_mRNA"), i, j] = (
                    sol.u[j][V.cfosmRNAc]
                )
                simulations[observables_index("cFos_Protein"), i, j] = (
                    (sol.u[j][V.pcFOSn] + sol.u[j][V.cFOSn]) * (p[C.Vn] / p[C.Vc])
                    + sol.u[j][V.cFOSc] + sol.u[j][V.pcFOSc]
                )
                simulations[observables_index("Phosphorylated_cFos"), i, j] = (
                    sol.u[j][V.pcFOSn] * (p[C.Vn] / p[C.Vc]) + sol.u[j][V.pcFOSc]
                )
            end
        end
    end
end
end # module
