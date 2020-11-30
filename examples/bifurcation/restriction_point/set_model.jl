function diffeq(du,u,p,t)
    du[V.M] = p[C.kM]*p[C.S]/(p[C.KS]+p[C.S]) - p[C.dM]*u[V.M]
    du[V.E] = p[C.kE]*(u[V.M]/(p[C.KM]+u[V.M]))*(u[V.E]/(p[C.KE]+u[V.E])) + p[C.kb]*u[V.M]/(p[C.KM]+u[V.M]) +
                p[C.kP1]*u[V.CD]*u[V.RE]/(p[C.KCD]+u[V.RE]) + p[C.kP2]*u[V.CE]*u[V.RE]/(p[C.KCE]+u[V.RE]) - p[C.dE]*u[V.E]-p[C.kRE]*u[V.R]*u[V.E]
    du[V.CD] = p[C.kCD]*u[V.M]/(p[C.KM]+u[V.M]) + p[C.kCDS]*p[C.S]/(p[C.KS]+p[C.S]) - p[C.dCD]*u[V.CD]
    du[V.CE] = p[C.kCE]*u[V.E]/(p[C.KE]+u[V.E]) - p[C.dCE]*u[V.CE]
    du[V.R] = p[C.kR] + p[C.kDP]*u[V.RP]/(p[C.KRP]+u[V.RP]) - p[C.kRE]*u[V.R]*u[V.E] - p[C.kP1]*u[V.CD]*u[V.R]/(p[C.KCD]+u[V.R]) -
                p[C.kP2]*u[V.CE]*u[V.R]/(p[C.KCE]+u[V.R]) - p[C.dR]*u[V.R]
    du[V.RP] = p[C.kP1]*u[V.CD]*u[V.R]/(p[C.KCD]+u[V.R]) + p[C.kP2]*u[V.CE]*u[V.R]/(p[C.KCE]+u[V.R]) + p[C.kP1]*u[V.CD]*u[V.RE]/(p[C.KCD]+u[V.RE]) +
                p[C.kP2]*u[V.CE]*u[V.RE]/(p[C.KCE]+u[V.RE]) - p[C.kDP]*u[V.RP]/(p[C.KRP]+u[V.RP]) - p[C.dRP]*u[V.RP]
    du[V.RE] = p[C.kRE]*u[V.R]*u[V.E]-p[C.kP1]*u[V.CD]*u[V.RE]/(p[C.KCD]+u[V.RE]) - p[C.kP2]*u[V.CE]*u[V.RE]/(p[C.KCE]+u[V.RE]) - p[C.dRE]*u[V.RE]

end


function param_values()::Vector{Float64}
    p::Vector{Float64} = zeros(C.NUM)

    p[C.S] = 2.0
    p[C.kE] = 0.4
    p[C.kM] = 1.0
    p[C.kCD] = 0.03
    p[C.kCDS] = 0.45
    p[C.kR] = 0.18
    p[C.kRE] = 180
    p[C.kb] = 0.003
    p[C.KS] = 0.5
    p[C.kCE] = 0.35
    p[C.dM] = 0.7
    p[C.dE] = 0.25
    p[C.dCD] = 1.5
    p[C.dCE] = 1.5
    p[C.dR] = 0.06
    p[C.dRP] = 0.06
    p[C.dRE] = 0.03
    p[C.kP1] = 18.0
    p[C.kP2] = 18.0
    p[C.kDP] = 3.6
    p[C.KM] = 0.15
    p[C.KE] = 0.15
    p[C.KCD] = 0.92
    p[C.KCE] = 0.92
    p[C.KRP] = 0.01

    return p
end


function get_derivatives(u::Vector{Float64},p::Vector{Float64})
    # derivatives: dF/d[bifurcation_param]
    dFdp::Vector{Float64} = zeros(V.NUM)

    dFdp[V.M] = p[C.kM]*p[C.KS]/(p[C.KS] + p[C.S])^2
    dFdp[V.CD] = p[C.kCDS]*p[C.KS]/(p[C.KS] + p[C.S])^2

    return dFdp
end


function get_steady_state(p::Vector{Float64})
    tspan::Tuple{Float64,Float64} = (0.0,Inf)
    u0::Vector{Float64} = zeros(V.NUM)

    prob = ODEProblem(diffeq,u0,tspan,p)
    prob = SteadyStateProblem(prob)
    sol = solve(prob,DynamicSS(CVODE_BDF()),dt=1.0)

    return sol.u
end