@testset "Parameter Estimation" begin
    let model_ode = Model("../examples/fos_model")
        output = []
        @testset "optimization" begin
            initpop = generate_initial_population(model_ode)
            @test size(initpop) == (225, 75)
            scipy_differential_evolution(model_ode, 1, maxiter=10, init=initpop)
            lines = open(joinpath(model_ode.path, "fitparam", "1", "optimization.log"), "r") do f
                readlines(f)
            end
            @test startswith(lines[end], "differential_evolution step 10:")
            push!(output, "logs")
            push!(output, "fitparam")
        end
        @testset "visualization" begin
            @test run_simulation(model_ode, viz_type="best") === nothing
            files = readdir(joinpath(model_ode.path, "figure", "simulation", "best"))
            n_pdf = 0
            for file in files
                if occursin(".pdf", file)
                    n_pdf += 1
                end
            end
            @test n_pdf == length(model_ode.observables)
            push!(output, "figure")
        end
        for dir in output
            rm(joinpath(model_ode.path, "$dir"), recursive=true, force=true)
        end
    end
end
