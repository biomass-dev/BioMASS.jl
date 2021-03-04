function diffeq!(du,u,p,t)
    CycE = p[C.CycElevel] - u[V.CycEp27]

    Vdp27 = p[C.kd27] + (p[C.kd27e]*CycE)*u[V.Skp2]
    Vdcyce = p[C.kdcyce] + p[C.kdcycee]*CycE
    Vdskp2 = p[C.kdskp2] + p[C.kdskp2c1]*u[V.Cdh1]

    Vicdh1 = p[C.kicdh1e]*CycE

    du[V.p27T] = p[C.ks27] - Vdp27*u[V.p27T]

    du[V.Skp2] = p[C.ksskp2] - Vdskp2*u[V.Skp2]

    du[V.CycEp27] = p[C.kasse]*(p[C.CycElevel]-u[V.CycEp27])*(u[V.p27T]-u[V.CycEp27])-(p[C.kdise]+Vdp27+Vdcyce)*u[V.CycEp27]

    du[V.EmiC] = p[C.kasec]*(p[C.Cdh1T]-u[V.EmiC])*(p[C.Emi1T]-u[V.EmiC]) - (p[C.kdiec]+p[C.kdemi1])*u[V.EmiC]

    du[V.Cdh1dp] = p[C.kacdh1]*(p[C.Cdh1T]-u[V.Cdh1dp]) - Vicdh1*u[V.Cdh1dp]

    du[V.Cdh1] = (p[C.kdiec]+p[C.kdemi1])*(u[V.Cdh1dp]-u[V.Cdh1]) - p[C.kasec]*u[V.Cdh1]*(p[C.Emi1T]-u[V.EmiC])+p[C.kacdh1]*(p[C.Cdh1T]-u[V.EmiC]-u[V.Cdh1])-Vicdh1*u[V.Cdh1]

end


function param_values()::Vector{Float64}
    p::Vector{Float64} = zeros(C.NUM)

    p[C.kscyce] = 0.003
    p[C.kdcyce] = 0.001
    p[C.kdcycee] = 0.0001
    p[C.kdcycea] = 0.03
    p[C.kasse] = 1
    p[C.kdise] = 0.02
    ## CYCA SYNTHESISp[C.DEGRADATION AND P27 BINDING/DISSOCIATION:
    p[C.kscyca] = 0.0025
    p[C.kdcyca] = 0.002
    p[C.kdcycac1] = 0.4
    p[C.kassa] = 1
    p[C.kdisa] = 0.02
    ## P27 SYNTHESIS AND DEGRADATION:
    p[C.ks27] = 0.008
    p[C.kd27] = 0.004
    p[C.kd27e] = 2
    p[C.kd27a] = 2
    ## EMI1 SYNTHESIS AND DEGRADATION:
    p[C.ksemi1] = 0.003
    p[C.kdemi1] = 0.001
    ## CDH1 REGULATION:
    p[C.Cdh1T] = 1
    p[C.kacdh1] = 0.02
    p[C.kicdh1e] = 0.07
    p[C.kicdh1a] = 0.2
    p[C.kasec] = 2
    p[C.kdiec] = 0.02
    ## SKP2 SYNTHESIS AND DEGRADATION:
    p[C.ksskp2] = 0.004
    p[C.kdskp2] = 0.002
    p[C.kdskp2c1] = 0.2
    ## CDK INHIBITOR
    p[C.Inhibitor] = 0.0

    p[C.Emi1T] = 0.0
    p[C.CycElevel] = 1.0

    return p
end


function get_derivatives(u::Vector{Float64},p::Vector{Float64})
    # derivatives: dF/d[bifurcation_param]
    dFdp::Vector{Float64} = zeros(V.NUM)

    dFdp[V.p27T]  = -p[C.kd27e]*u[V.Skp2]*u[V.p27T]
    dFdp[V.CycEp27] = p[C.kasse]*(1.0-u[V.CycEp27])*(u[V.p27T]-u[V.CycEp27])-(p[C.kd27e]*u[V.Skp2]+p[C.kdcycee])*u[V.CycEp27]
    dFdp[V.Cdh1dp] = -p[C.kicdh1e]*u[V.Cdh1dp]
    dFdp[V.Cdh1] =  -p[C.kicdh1e]*u[V.Cdh1]

    return dFdp
end


function get_steady_state(p::Vector{Float64})
    tspan::Tuple{Float64,Float64} = (0.0,Inf)
    u0::Vector{Float64} = zeros(V.NUM)

    prob = ODEProblem(diffeq!,u0,tspan,p)
    prob = SteadyStateProblem(prob)
    sol = solve(prob,DynamicSS(CVODE_BDF()),dt=1.0)

    return sol.u
end