using Documenter, SchoofPointCounts

DocMeta.setdocmeta!(SchoofPointCounts, :DocTestSetup, :(using SchoofPointCounts); recursive=true)

makedocs(;
    modules=[SchoofPointCounts],
    authors="William Herrera",
    sitename="SchoofPointCounts.jl Documentation",
    repo=Documenter.Remotes.GitHub("wherrera10", "SchoofPointCounts.jl"),
    format=Documenter.HTML(;
        canonical="https://wherrera10.github.io/SchoofPointCounts.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        # Add future pages here, e.g., "API Reference" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/wherrera10/SchoofPointCounts.jl.git",
    devbranch="main",
)

  
