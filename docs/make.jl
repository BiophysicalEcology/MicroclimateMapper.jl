using Documenter, DocumenterVitepress
using MicroclimateMapper

makedocs(;
    modules=[MicroclimateMapper],
    authors="Michael Kearney, Rafael Schouten",
    sitename="MicroclimateMapper.jl",
    clean=true,
    doctest=true,
    checkdocs=:none,
    format=DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/BiophysicalEcology/MicroclimateMapper.jl",
        devbranch = "main",
        devurl = "dev";
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "Raster" => "examples/raster.md",
            "Vector" => "examples/vector.md",
        ],
        "API" => "api.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo="github.com/BiophysicalEcology/MicroclimateMapper.jl",
    branch="gh-pages",
    target = joinpath(@__DIR__, "build"),
    devbranch="main",
    push_preview=true,
)
