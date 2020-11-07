using BioMASS
using Test

@time begin
    @testset "BioMASS.jl" begin
        include("parameter_estimation.jl")
    end
end
