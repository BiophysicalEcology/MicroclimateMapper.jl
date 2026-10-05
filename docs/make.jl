using Documenter
using DocumenterVitepress
using MicroclimateMapper
using CairoMakie

# Don't output huge svgs for Makie plots
CairoMakie.activate!(type = "png")

# Downloaded data sets go here unless a location is already set
haskey(ENV, "RASTERDATASOURCES_PATH") || (ENV["RASTERDATASOURCES_PATH"] = mkpath(joinpath(@__DIR__, "data_cache")))

# Helpers for the figures, loaded in the examples with `using Main.FigureHelpers`
include("figure_helpers.jl")

makedocs(
    modules = [MicroclimateMapper],
    sitename = "MicroclimateMapper.jl",
    authors = "Michael Kearney, Rafael Schouten et al.",
    clean = true,
    doctest = false,
    checkdocs = :exports,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/BiophysicalEcology/MicroclimateMapper.jl", # this must be the full URL!
        devbranch = "main",
        devurl = "dev";
    ),
    source = "src",
    build = "build",
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/BiophysicalEcology/MicroclimateMapper.jl",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
