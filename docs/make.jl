using BioMASS
using Documenter

makedocs(;
    modules=[BioMASS],
    authors="Hiroaki Imoto <himoto@protein.osaka-u.ac.jp>",
    repo="https://github.com/himoto/BioMASS.jl/blob/{commit}{path}#L{line}",
    sitename="BioMASS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://himoto.github.io/BioMASS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/himoto/BioMASS.jl",
)
