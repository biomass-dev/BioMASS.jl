import BioMASS:isinstalled

@testset "Parameter Estimation" begin
    model_ode = load_model("../examples/fos_model")
    output = []
    @testset "optimization" begin
        optimize(
            model_ode, 1, max_generation=3, popsize=3,
            local_search_method="mutation", n_children=15
        )
        lines = open(joinpath(model_ode.path, "logs", "1.log"), "r") do f
            readlines(f)
        end
        @test lines[end][1:13] == "Generation3: "
        push!(output, "logs")
        push!(output, "fitparam")
    end
    @testset "optimization_continue" begin
        optimize_continue(
            model_ode, 1, max_generation=6, popsize=3,
            local_search_method="CMAES", maxiter=30
        )
        lines = open(joinpath(model_ode.path, "logs", "1.log"), "r") do f
            readlines(f)
        end
        @test lines[end][1:13] == "Generation6: "
        # test differential_evolution
        if isinstalled("scipy.optimize")
            @testset "Differential evolution" begin
                optimize_continue(
                    model_ode, 1, max_generation=9, popsize=3,
                    local_search_method="DE", maxiter=10
                )
                lines = open(joinpath(model_ode.path, "logs", "1.log"), "r") do f
                    readlines(f)
                end
                @test lines[end][1:13] == "Generation9: "
            end
            @testset "Modified Powell's method" begin
                optimize_continue(
                    model_ode, 1, max_generation=12, popsize=3,
                    local_search_method="Powell", maxiter=5
                )
                lines = open(joinpath(model_ode.path, "logs", "1.log"), "r") do f
                    readlines(f)
                end
                @test lines[end][1:14] == "Generation12: "
            end
        end
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