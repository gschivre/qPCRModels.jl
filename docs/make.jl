using qPCRModels
using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

DocMeta.setdocmeta!(qPCRModels, :DocTestSetup, :(using qPCRModels); recursive = true)

makedocs(;
    plugins = [bib],
    checkdocs = :exports,
    modules = [qPCRModels],
    authors = "Geoffrey Schivre",
    sitename = "qPCRModels.jl",
    format = Documenter.HTML(;
        canonical = "https://gschivre.github.io/qPCRModels.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "The Birch model" => "birch.md",
        "API reference" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/gschivre/qPCRModels.jl",
    devbranch = "main",
)
