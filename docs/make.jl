using qPCRModels
using Documenter

# The home and Birch model pages are generated from README.md so the two never drift apart.
# README math uses $$...$$ (display, on their own lines) and $...$ (inline), which GitHub renders;
# Documenter wants ```math blocks and ``...`` inline math.
function readme_section(readme, from, to)
    i = findfirst(from, readme)
    j = findnext(to, readme, last(i))
    return readme[first(i):prevind(readme, first(j))]
end
function to_documenter(md)
    md = replace(md, r"^([ \t]*)\$\$\n(.*?)\n[ \t]*\$\$"ms => s"\1```math\n\2\n\1```")
    md = replace(md, r"\$(.+?)\$" => s"``\1``")
    md = replace(md, "](docs/src/assets/" => "](assets/") # generated pages live in docs/src
    return replace(md, "](#discrete-cycles-vs-continuous-time)" => "](@ref \"Discrete cycles vs continuous time\")")
end

let readme = read(joinpath(@__DIR__, "..", "README.md"), String)
    intro = readme_section(readme, "`qPCRModels.jl` is a Julia package", "## Why the Birch model?")
    usage = readme_section(readme, "## Installation", "## Package layout")
    birch = readme_section(readme, "## Why the Birch model?", "## Installation")
    write(joinpath(@__DIR__, "src", "index.md"), to_documenter("""
        ```@meta
        CurrentModule = qPCRModels
        ```

        # qPCRModels.jl

        $(intro)See [The Birch model](@ref) for the mechanistic derivation of the default model.

        $(usage)"""))
    write(joinpath(@__DIR__, "src", "birch.md"), to_documenter("# The Birch model\n\n" * birch))
end

DocMeta.setdocmeta!(qPCRModels, :DocTestSetup, :(using qPCRModels); recursive=true)

makedocs(;
    checkdocs=:exports,
    modules=[qPCRModels],
    authors="Geoffrey Schivre",
    sitename="qPCRModels.jl",
    format=Documenter.HTML(;
        canonical="https://gschivre.github.io/qPCRModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "The Birch model" => "birch.md",
        "API reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/gschivre/qPCRModels.jl",
    devbranch="main",
)
