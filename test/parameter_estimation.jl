import BioMASS: isinstalled

const model = load_model("../examples/fos_model")

output = []

@testset "Parameter Estimation" begin
    @testset "optimization" begin
        optimize(
            model, 1, max_generation=3, popsize=3,
            local_search_method="mutation", n_children=15
        )
        lines = open(model.path * "/logs/1.log", "r") do f
            readlines(f)
        end
        @test lines[end][1:13] == "Generation3: "
        push!(output, "logs")
        push!(output, "fitparam")
    end
    @testset "optimization_continue" begin
        if isinstalled("scipy.optimize")
            optimize_continue(
                model, 1, max_generation=6, popsize=3,
                local_search_method="Powell"
            )
        else
            optimize_continue(
                model, 1, max_generation=6, popsize=3,
                local_search_method="mutation", n_children=15
            )
        end
        lines = open(model.path * "/logs/1.log", "r") do f
            readlines(f)
        end
        @test lines[end][1:13] == "Generation6: "
        # test differential_evolution
        if isinstalled("scipy.optimize")
            optimize_continue(
                model, 1, max_generation=9, popsize=3,
                local_search_method="DE"
            )
            lines = open(model.path * "/logs/1.log", "r") do f
                readlines(f)
            end
            @test lines[end][1:13] == "Generation9: "
        end
    end
    if isinstalled("matplotlib")
        @testset "visualization" begin
            visualize(model, viz_type="best")
            files = readdir("../examples/fos_model/figure/simulation/best")
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
            @test param2biomass(model.path) === nothing
            push!(output, "dat2npy")
        end
    end
    for dir in output
        rm("../examples/fos_model/$dir", recursive=true, force=true)
    end
end