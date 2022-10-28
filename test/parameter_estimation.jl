import BioMASS: isinstalled
using PyCall

@testset "Parameter Estimation" begin
    model_ode = Model("../examples/fos_model")
    output = []
    @testset "optimization" begin
        # initpop = generate_initial_population(model_ode);
        optimizer_options = Dict{String, Int}(py"{'maxiter': 10}");
        scipy_differential_evolution(model_ode, 1, optimizer_options);
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
    if isinstalled("numpy")
        @testset "conversion" begin
            @test param2biomass(model_ode.path) === nothing
            @test isdir(joinpath(model_ode.path, "dat2npy", "out", "1"))
            push!(output, "dat2npy")
        end
    end
    for dir in output
        rm(joinpath(model_ode.path, "$dir"), recursive=true, force=true)
    end
end
