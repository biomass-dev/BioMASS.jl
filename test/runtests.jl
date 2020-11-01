using BioMASS
using Test

const MODEL_PATH = "../fos_model"
include(MODEL_PATH * "/exec_model.jl")

using .ExecModel

@testset "BioMASS.jl" begin
    optimize(MODEL_PATH, 1, max_generation=10)
    optimize_continue(MODEL_PATH, 1, max_generation=20)
    # param2biomass(MODEL_PATH)
end
