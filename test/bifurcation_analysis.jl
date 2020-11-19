using DelimitedFiles
using Sundials
using SteadyStateDiffEq
import BioMASS: isinstalled

@testset "Bifurcation analysis" begin
    for model in ["restriction_point", "g1s_transition", "mitotic_exit"]
        @testset "$model" begin
            MODEL_PATH = "../examples/bifurcation/" * model
            create_diffeq(MODEL_PATH)
            include(MODEL_PATH * "/diagram.jl")
            calc_fixed_point_vec(MODEL_PATH)
            for file in ["/data/fp.dat", "/data/ev.dat"]
                @test isfile(MODEL_PATH * file)
            end
            if isinstalled("matplotlib")
                using PyPlot
                bifurcation_diagram(MODEL_PATH)
                @test isfile(MODEL_PATH * "/bifurcation_diagram.pdf")
            end
            rm(MODEL_PATH * "/forwarddiff.jl")
            rm(MODEL_PATH * "/data", recursive=true, force=true)
        end
    end
end