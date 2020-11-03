using BioMASS
using Test

const MODEL_PATH = "../fos_model"
include(MODEL_PATH * "/exec_model.jl")

using .ExecModel

@testset "BioMASS.jl" begin
    optimize(MODEL_PATH, 1, max_generation=10)
    lines = open(MODEL_PATH * "/logs/1.log", "r") do f
        readlines(f)
    end
    @test lines[end][1:14] == "Generation10: "

    optimize_continue(MODEL_PATH, 1, max_generation=20)
    lines = open(MODEL_PATH * "/logs/1.log", "r") do f
        readlines(f)
    end
    @test lines[end][1:14] == "Generation20: "
    
    param2biomass(MODEL_PATH)
end
