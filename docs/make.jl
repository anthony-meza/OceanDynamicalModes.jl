using OceanDynamicalModes
using Documenter

DocMeta.setdocmeta!(OceanDynamicalModes, :DocTestSetup, :(using OceanDynamicalModes); recursive=true)

makedocs(;
    modules=[OceanDynamicalModes],
    authors="Anthony Meza",
    repo="https://github.com/anthony-meza/OceanDynamicalModes.jl/blob/{commit}{path}#{line}",
    sitename="OceanDynamicalModes.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://anthony-meza.github.io/OceanDynamicalModes.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/anthony-meza/OceanDynamicalModes.jl",
    devbranch="main",
)
