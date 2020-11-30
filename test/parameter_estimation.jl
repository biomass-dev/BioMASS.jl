import BioMASS: isinstalled

const model = load_model("../examples/fos_model")

output = []

@testset "Parameter Estimation" begin
    @testset "optimization" begin
        optimize(model, 1, max_generation=10)
        lines = open(model.path * "/logs/1.log", "r") do f
            readlines(f)
        end
        @test lines[end][1:14] == "Generation10: "
        push!(output, "logs")
        push!(output, "fitparam")
    end
    @testset "optimization_continue" begin
        optimize_continue(model, 1, max_generation=20)
        lines = open(model.path * "/logs/1.log", "r") do f
            readlines(f)
        end
        @test lines[end][1:14] == "Generation20: "
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