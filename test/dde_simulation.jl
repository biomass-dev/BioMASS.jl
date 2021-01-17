import BioMASS: isinstalled


@testset "DDE simulation" begin
    model_dde = load_model("../examples/nfkb_model")
    if isinstalled("matplotlib")
        @test visualize(model_dde, viz_type="original") === nothing
        rm("../examples/nfkb_model/figure", recursive=true, force=true)
    end
end