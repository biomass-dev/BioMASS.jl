import BioMASS: isinstalled
using PyCall

@testset "Parameter Estimation" begin
    model_ode = Model("../examples/fos_model")
    output = []
    @testset "optimization" begin
        initpop = generate_initial_population(model_ode)
        scipy_differential_evolution(model_ode, 1, maxiter=10, init=initpop)
        lines = open(joinpath(model_ode.path, "fitparam", "1", "optimization.log"), "r") do f
            readlines(f)
        end
        @test startswith(lines[end], "differential_evolution step 10:")
        push!(output, "logs")
        push!(output, "fitparam")
    end
    if isinstalled("matplotlib")
        @testset "visualization" begin
            @test run_simulation(model_ode, viz_type="best") === nothing
            files = readdir(joinpath(model_ode.path, "figure", "simulation", "best"))
            n_pdf = 0
            for file in files
                if occursin(".pdf", file)
                    n_pdf += 1
                end
            end
            @test n_pdf == 8  # length(observables)
            push!(output, "figure")
        end
    end
    for dir in output
        rm(joinpath(model_ode.path, "$dir"), recursive=true, force=true)
    end
end
