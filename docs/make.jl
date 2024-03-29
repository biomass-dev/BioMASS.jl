using BioMASS
using Documenter

makedocs(;
    modules=[BioMASS],
    authors="Hiroaki Imoto <hiroaki.imoto@ucd.ie>",
    repo="https://github.com/biomass-dev/BioMASS.jl/blob/{commit}{path}#L{line}",
    sitename="BioMASS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://biomass-dev.github.io/BioMASS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting started with BioMASS.jl" => [
            "Parameter Estimation" => "usage/parameter_estimation.md",
            "Bifurcation Analysis" => "usage/bifurcation_analysis.md",
        ],
        "References" => "references.md"
    ],
)

deploydocs(;
    repo="github.com/biomass-dev/BioMASS.jl",
)
