using BioMASS
using Test

@time begin
    @testset "BioMASS.jl" begin
        include("parameter_estimation.jl")
        # include("dde_simulation.jl")
        include("bifurcation_analysis.jl")
    end
end
