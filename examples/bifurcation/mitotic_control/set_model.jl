function diffeq!(du,u,p,t)
    Complex = p[C.B55tot] - u[V.PP2AB55]
    Wee1p = 1 - u[V.Wee1] - u[V.Wee1pp]
    Cdc25p = 1 - u[V.Cdc25] - u[V.Cdc25pp]
    Vwee = (p[C.kweeS]*(1-u[V.Wee1]) + p[C.kweeF]*u[V.Wee1])
    V25 = p[C.k25S]*(1-u[V.Cdc25pp]) + p[C.k25F]*u[V.Cdc25pp]
    VGwl = p[C.kGwENSA]*u[V.Gwlp]

    du[V.Subp] = p[C.kcBc1Sub]*u[V.CycBCdk1]/((1 + (p[C.InhCDK]/p[C.Kd])))*(p[C.SubT]-u[V.Subp]) - p[C.kB55Sub]*u[V.PP2AB55]*u[V.Subp]
    du[V.CycBCdk1] = V25*(p[C.CycBCdk1T] - u[V.CycBCdk1]) - Vwee*u[V.CycBCdk1]
    du[V.PP1] = (p[C.kapp1] + p[C.kapp1a]*u[V.PP1])*(p[C.PP1T] - u[V.PP1]) - (p[C.kipp1] + p[C.kipp1C]*u[V.CycBCdk1]/((1 + (p[C.InhCDK]/p[C.Kd]))))*u[V.PP1]
    du[V.pENSAt] = VGwl*(p[C.ENSAtot] - u[V.pENSAt]) - p[C.kcatB55]*Complex
    du[V.Gwlp] = (p[C.kcBc1G]*u[V.CycBCdk1]/((1 + (p[C.InhCDK]/p[C.Kd]))) + p[C.kcAc2G]*p[C.CycACdk2T])*(p[C.Gwtot] - u[V.Gwlp]) - (p[C.kB55G]*u[V.PP2AB55] + p[C.kppxGwl] + p[C.kPP1Gw]*u[V.PP1])*u[V.Gwlp]
    du[V.PP2AB55] = p[C.kdis]*Complex + p[C.kcatB55]*Complex - p[C.kass]*u[V.PP2AB55]*(u[V.pENSAt] - Complex)
    du[V.Wee1] = (p[C.kppxY15] + p[C.kB55W1]*u[V.PP2AB55])*Wee1p - (p[C.kcBc1W1]*u[V.CycBCdk1]/((1 + (p[C.InhCDK]/p[C.Kd]))) + p[C.kcAc2W1]*p[C.CycACdk2T])*u[V.Wee1]
    du[V.Wee1pp] = (p[C.kcBc1W1]*u[V.CycBCdk1]/((1 + (p[C.InhCDK]/p[C.Kd]))) + p[C.kcAc2W1]*p[C.CycACdk2T])*Wee1p - (p[C.kppxY15] + p[C.kB55W1]*u[V.PP2AB55])*u[V.Wee1pp]
    du[V.Cdc25] = (p[C.kppxY15] + p[C.kB5525]*u[V.PP2AB55])*Cdc25p - (p[C.kcBc125]*u[V.CycBCdk1]/((1 + (p[C.InhCDK]/p[C.Kd]))) + p[C.kcAc225]*p[C.CycACdk2T])*u[V.Cdc25]
    du[V.Cdc25pp] = (p[C.kcBc125]*u[V.CycBCdk1]/((1 + (p[C.InhCDK]/p[C.Kd]))) + p[C.kcAc225]*p[C.CycACdk2T])*(Cdc25p) - (p[C.kppxY15] + p[C.kB5525]*u[V.PP2AB55])*u[V.Cdc25pp]

end


function param_values()::Vector{Float64}
    p::Vector{Float64} = zeros(C.NUM)

    p[C.InhCDK] = 2.0
    p[C.CycBCdk1T] = 8.1808
    p[C.CycACdk2T] = 1.0000
    p[C.PP1T] = 1.0000
    p[C.kapp1] = 0.0115
    p[C.kapp1a] = 0.7054
    p[C.kipp1] = 0.0018
    p[C.kipp1C] = 0.7549
    p[C.kPP1Gw] = 18.4724
    p[C.ENSAtot] = 1.0000
    p[C.B55tot] = 0.2500
    p[C.SubT] = 1.0000
    p[C.kass] = 617.2807
    p[C.kdis] = 0.0088
    p[C.kcatB55] = 1.0338
    p[C.kGwENSA] = 20.8811
    p[C.kppxGwl] = 0.1560
    p[C.kcBc1Sub] = 0.0080
    p[C.kcBc1G] = 0.2393
    p[C.Gwtot] = 1.0000
    p[C.kB55G] = 496.5636
    p[C.kB55Sub] = 0.0593
    p[C.kcAc2G] = 0.1916
    p[C.k25S] = 0.0050
    p[C.k25F] = 0.9411
    p[C.kweeS] = 0.0050
    p[C.kweeF] = 47.2937
    p[C.kcBc1W1] = 1.3132
    p[C.kcBc125] = 1.3132
    p[C.kppxY15] = 0.0050
    p[C.kcAc2W1] = 0.1096
    p[C.kcAc225] = 0.1096
    p[C.kB55W1] = 0.5511
    p[C.kB5525] = 0.5511
    p[C.Kd] = 0.025

    return p
end


function get_derivatives(u::Vector{Float64},p::Vector{Float64})
    # derivatives: dF/d[bifurcation_param]
    dFdp::Vector{Float64} = zeros(V.NUM)

    Wee1p = 1 - u[V.Wee1] - u[V.Wee1pp]
    Cdc25p = 1 - u[V.Cdc25] - u[V.Cdc25pp]

    dFdp[V.Subp] = - (1/p[C.Kd])*(p[C.kcBc1Sub]*u[V.CycBCdk1]/(1 + (p[C.InhCDK]/p[C.Kd]))^2)*(p[C.SubT]-u[V.Subp])
    dFdp[V.PP1] = (1/p[C.Kd])*(p[C.kipp1] + p[C.kipp1C]*u[V.CycBCdk1]/(1 + (p[C.InhCDK]/p[C.Kd]))^2)*u[V.PP1]
    dFdp[V.Gwlp] = - (1/p[C.Kd])*(p[C.kcBc1G]*u[V.CycBCdk1]/(1 + (p[C.InhCDK]/p[C.Kd]))^2)*(p[C.Gwtot] - u[V.Gwlp])
    dFdp[V.Wee1] = (1/p[C.Kd])*(p[C.kcBc1W1]*u[V.CycBCdk1]/(1 + (p[C.InhCDK]/p[C.Kd]))^2)*u[V.Wee1]
    dFdp[V.Wee1pp] = - (1/p[C.Kd])*(p[C.kcBc1W1]*u[V.CycBCdk1]/(1 + (p[C.InhCDK]/p[C.Kd]))^2)*Wee1p
    dFdp[V.Cdc25] = (1/p[C.Kd])*(p[C.kcBc125]*u[V.CycBCdk1]/(1 + (p[C.InhCDK]/p[C.Kd]))^2)*u[V.Cdc25]
    dFdp[V.Cdc25pp] = - (1/p[C.Kd])*(p[C.kcBc125]*u[V.CycBCdk1]/(1 + (p[C.InhCDK]/p[C.Kd]))^2)*(Cdc25p)

    return dFdp
end


function get_steady_state(p::Vector{Float64})
    tspan::Tuple{Float64,Float64} = (0.0,Inf)
    u0::Vector{Float64} = zeros(V.NUM)
    u0[V.PP1] = 1.0
    u0[V.PP2AB55] = 0.25
    u0[V.Wee1] = 1.0
    u0[V.Cdc25] = 1.0

    prob = ODEProblem(diffeq!,u0,tspan,p)
    prob = SteadyStateProblem(prob)
    sol = solve(prob,DynamicSS(CVODE_BDF()),dt=1.0)

    return sol.u
end