const model = load_model("../fos_model")

@testset "Parameter Estimation" begin
    @testset "optimize" begin
        optimize(model, 1, max_generation=10)
        lines = open(model.path * "/logs/1.log", "r") do f
            readlines(f)
        end
        @test lines[end][1:14] == "Generation10: "
    end
    @testset "optimize_continue" begin
        optimize_continue(model, 1, max_generation=20)
        lines = open(model.path * "/logs/1.log", "r") do f
            readlines(f)
        end
        @test lines[end][1:14] == "Generation20: "
    end
    @testset "conversion" begin
        @test param2biomass(model.path) === nothing
    end
end