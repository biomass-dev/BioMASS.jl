using ForwardDiff:jacobian

const MC = 100000           # maximum of counts
const IVAL = 1e-2           # first variation
const RATE = 1e-3           # variation rate
const NEPS = 1e-12          # eps of Newton's method


function create_diffeq(model_path::String)
    lines::Vector{String} = []
    open(joinpath(model_path, "set_model.jl"), "r") do f
        append!(lines, readlines(f))
    end
    for (i, line) in enumerate(lines)
        if occursin("function diffeq!", line)
            lines[i] = "function diffeq(u::Vector)\n    du = similar(u)\n\n"
        elseif line == "end"
            lines[i] = "    return du\nend"
            lines = lines[1:i]
            break
        end
    end
    open(joinpath(model_path, "forwarddiff.jl"), "w") do f
        for line in lines
            write(f, line * "\n")
        end
    end
end


# matrix transformation (large diagonal elements move to upper)
function pivoting!(s::Matrix{Float64}, pivot::Int, dim_newton::Int)
    v0::Vector{Float64} = zeros(dim_newton + 1)
    v1::Vector{Float64} = zeros(dim_newton + 1)
    possess::Int = 0
    max_element::Float64 = 0.0

    for i in pivot:size(s, 1)
        current_element = abs(s[i, pivot])
        if max_element <= current_element
            max_element = current_element
            possess = i
        end
    end

    for j in 1:size(s, 2)
        v0[j] = s[possess, j]
        v1[j] = s[pivot, j]
    end

    for j in 1:size(s, 2)
        s[possess, j] = v1[j]
        s[pivot, j] = v0[j]
    end
end

# Gaussian elimination (row reduction)
function gaussian_elimination!(s::Matrix{Float64}, e::Vector{Float64}, dim_newton::Int)
    for i in 1:dim_newton
        pivoting!(s, i, dim_newton)
    end
    # forward
    for k in 1:size(s, 1)
        w = (s[k, k] != 0.0) ? 1.0 / s[k, k] : 1.0
        for j in k:size(s, 2)
            s[k, j] *= w
            for i in k:size(s, 1)
                s[i, j] -= s[i, k] * s[k, j]
            end
        end
    end
    # backward
    for i in size(s, 1):-1:1
        sum = 0.0
        for j in i:length(e)
            sum += s[i, j] * e[j]
        end
        e[i] = s[i, end] - sum
    end
end


# Newton's method
function newtons_method!(
        diffeq::Function,
        get_derivatives::Function,
        x::Vector{Float64},
        real_part::Vector{Float64},
        imaginary_part::Vector{Float64}, 
        fix_num::Int,
        p::Vector{Float64},
        successful::Bool,
        bifparam::Int,
        n_state::Int,
        dim_newton::Int,
        n_variable::Int)
    u::Vector{Float64} = zeros(n_state)
    vx::Vector{Float64} = zeros(dim_newton)
    s::Matrix{Float64} = zeros(dim_newton, dim_newton + 1)

    for i in eachindex(x)
        if fix_num == i
            for j in eachindex(vx)
                idx = i + j
                if idx > length(x)
                    idx -= length(x)
                end
                vx[j] = x[idx]
            end
            break
        else
            continue
        end
    end

    # initial error
    e::Vector{Float64} = zeros(dim_newton)
    error::Float64 = 1.0

    while error > NEPS
        for i in 1:n_variable
            if fix_num == i
                idx_param = n_variable - i
                p[bifparam] = (idx_param == 0) ? x[fix_num] : vx[idx_param]
                for j in eachindex(u)
                    idx = j - i
                    if idx == 0
                        u[j] = x[fix_num]
                    elseif idx < 0
                        u[j] = vx[n_variable + idx]
                    else
                        u[j] = vx[idx]
                    end
                end
                break
        else
                continue
        end
        end

        # initialization
        dFdx::Matrix{Float64} = jacobian(diffeq, u)
        dFdp::Vector{Float64} = get_derivatives(u, p)
        
        F::Vector{Float64} = diffeq(u)

        eigenvalues::Array{Complex{Float64},1} = eigvals(dFdx)
        for (i, eigenvalue) in enumerate(eigenvalues)
            real_part[i] = real(eigenvalue)
            imaginary_part[i] = imag(eigenvalue)
        end

        # s = [dF-F]
        for i in 1:n_variable
            if fix_num == i
                for k in 1:n_state
                    for j in 1:n_state
                        idx = i + j
                        if idx == n_variable
                            s[k, j] = dFdp[k]
                        elseif idx > n_variable
                            s[k, j] = dFdx[k, idx - n_variable]
                        else
                            s[k, j] = dFdx[k, idx]
                        end
                    end
                    s[k, n_variable] = -F[k]
                end
                break
            else
                continue
            end
        end

        gaussian_elimination!(s, e, dim_newton)

        # update error
        error = 0.0
        @inbounds for i in eachindex(e)
            vx[i] += e[i]
            error += e[i] * e[i]
        end
        error = sqrt(error)
        if isnan(error) || isinf(error)
            successful = false
            break
        end
    end

    for i in eachindex(x)
        if fix_num == i
            for j in eachindex(vx)
                idx = i + j
                if idx > length(x)
                    idx -= length(x)
                end
                x[idx] = vx[j]
            end
            break
        else
            continue
        end
    end
end


function new_curve!(
        model_path::Union{String,SubString{String}},
        p::Vector{Float64},
        diffeq::Function,
        get_derivatives::Function,
        get_steady_state::Function;
        direction::Bool=false,
        bifparam::Int,
        n_state::Int,
        n_param::Int=1,
        n_variable::Int=n_state + 1,
        dim_newton::Int=n_state)
    # bifparam : name(index) of bifurcation parameter
    # n_state : num of state variables
    # n_param : num of parameters
    # n_variable : num of variables
    # dim_newton : dim of Newton's method
    count::Int = 1
    x::Vector{Float64} = zeros(n_variable)
    dx::Vector{Float64} = zeros(n_variable)

    real_part::Vector{Float64} = zeros(n_state)
    imaginary_part::Vector{Float64} = zeros(n_state)

    # file
    if !isdir(
        joinpath(
            model_path,
            "data",
        )
    )
        mkdir(
            joinpath(
                model_path,
                "data",
            )
        )
    else
        files::Vector{String} = readdir(
            joinpath(
                model_path,
                "data",
            )
        )
        for file in files
            rm(
                joinpath(
                    model_path,
                    "data",
                    "$file",
                )
            )
        end
    end
    
    FOUT1 = open(joinpath(model_path, "data", "fp.dat"), "w") # file for fixed point
    FOUT2 = open(joinpath(model_path, "data", "ev.dat"), "w") # file for eigenvalues

    # initial condition
    x[1:n_state] = get_steady_state(p)
    x[end] = p[bifparam]  # x-axis

    # initial fixed
    fix_val::Float64 = x[end]
    fix_num::Int = n_variable
    x[fix_num] = fix_val

    # first Newton's method
    successful::Bool = true
    newtons_method!(
        diffeq, get_derivatives, x, real_part, imaginary_part, fix_num, p,
        successful, bifparam, n_state, dim_newton, n_variable
    )

    write(FOUT1, @sprintf("%d\t", count))
    for i in eachindex(x)
        write(FOUT1, @sprintf("%10.8e\t", x[i]))
    end
    write(FOUT1, @sprintf("%d\n", fix_num))
    write(FOUT2, @sprintf("%d\t", count))
    for i in 1:n_state
        write(
            FOUT2, @sprintf(
                "%10.8e\t%10.8e\t", real_part[i], imaginary_part[i]
            )
        )
    end
    write(FOUT2, @sprintf("%10.8e\t%d\n", p[bifparam], fix_num))
    count += 1

    # keep optimums
    px::Vector{Float64} = copy(x)

    # variation
    fix_val += ifelse(direction, +IVAL, -IVAL)

    # same fixed variable
    x[fix_num] = fix_val

    while count <= MC && successful
        newtons_method!(
            diffeq, get_derivatives, x, real_part, imaginary_part, fix_num, p,
            successful, bifparam, n_state, dim_newton, n_variable
        )

        # maximum variation
        for (i, prev) in enumerate(px)
        @inbounds dx[i] = x[i] - prev
        end
        sum::Float64 = 0.0
        for i in eachindex(dx)
            @inbounds sum += dx[i] * dx[i]
        end
        ave::Float64 = sqrt(sum)
        for i in eachindex(dx)
            @inbounds dx[i] /= ave
        end
        px = copy(x)
        for (i, diff) in enumerate(dx)
            @inbounds x[i] += abs(RATE) * diff
        end

        # fix variable with maximum variation
        fix_num = 1
        for i in 2:length(dx)
            if abs(dx[fix_num]) < abs(dx[i])
                fix_num = i
        end
        end

            # Stop calc.
        if x[end] <= 0.0
            successful = false
        end

        write(FOUT1, @sprintf("%d\t", count))
        for i in eachindex(x)
            write(FOUT1, @sprintf("%10.8e\t", x[i]))
        end
        write(FOUT1, @sprintf("%d\n", fix_num))
        write(FOUT2, @sprintf("%d\t", count))
        for i in 1:n_state
            write(
                FOUT2, @sprintf(
                    "%10.8e\t%10.8e\t", real_part[i], imaginary_part[i]
                )
            )
        end
        write(FOUT2, @sprintf("%10.8e\t%d\n", p[bifparam], fix_num))
        count += 1
    end

    close(FOUT1)
    close(FOUT2)
end


function get_bistable_regime(ev::Matrix{Float64}, n_state::Int)
    br::Vector{Int} = []
    for i in 1:size(ev, 1)
        if maximum(ev[i, [2j for j in 1:n_state]]) > 0.0
            push!(br, i)
        end
    end
    return br
end