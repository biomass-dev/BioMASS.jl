const MODEL_PATH = "../fos_model"
include(MODEL_PATH * "/exec_model.jl")

using .ExecModel

@testset "Parameter Estimation" begin
    @testset "optimize" begin
        optimize(MODEL_PATH, 1, max_generation=10)
        lines = open(MODEL_PATH * "/logs/1.log", "r") do f
            readlines(f)
        end
        @test lines[end][1:14] == "Generation10: "
    end

    @testset "optimize_continue" begin
        optimize_continue(MODEL_PATH, 1, max_generation=20)
        lines = open(MODEL_PATH * "/logs/1.log", "r") do f
            readlines(f)
        end
        @test lines[end][1:14] == "Generation20: "
    end

    @testset "conversion" begin
        @test param2biomass(MODEL_PATH) === nothing
    end
end