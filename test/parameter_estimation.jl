const model = load_model("../fos_model")

@testset "Parameter Estimation" begin
    @testset "optimization" begin
        optimize(model, 1, max_generation=10)
        lines = open(model.path * "/logs/1.log", "r") do f
            readlines(f)
        end
        @test lines[end][1:14] == "Generation10: "
    end
    @testset "visualization" begin
        visualize(model, viz_type="best")
        @test isdir("../fos_model/figure/simulation/best")
    end
    @testset "conversion" begin
        @test param2biomass(model.path) === nothing
    end
    for dir in ["logs", "fitparam", "figure", "fitparam"]
        rm("../fos_model/$dir", recursive=true, force=true)
    end
end