function diffeq(du,u,p,t)
    v = Dict{Int64,Float64}()
    

end


function param_values()::Vector{Float64}
    p::Vector{Float64} = zeros(C.NUM)

    return p
end


function initial_values()::Vector{Float64}
    u0::Vector{Float64} = zeros(V.NUM)

    return u0
end
