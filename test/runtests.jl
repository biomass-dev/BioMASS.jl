using BioMASS
using Test

@time begin
    @testset "BioMASS.jl" begin
        include("parameter_estimation.jl")
        include("bifurcation_analysis.jl")
    end
end
