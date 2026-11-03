# Build the standalone `qpcr_fit` executable from bin/qpcr_fit.jl with JuliaC (needs Julia ≥ 1.12):
#
#     julia --project=juliac -e 'using Pkg; Pkg.instantiate()'   # first time only
#     julia --project=juliac juliac/build.jl
#
# The bundle is written to juliac/qpcr_fit/: the executable is bin/qpcr_fit (bin/qpcr_fit.exe on Windows), and
# it needs the libraries next to it (lib/, or bin/ on Windows), but no Julia installation. The built-in datasets
# are copied to share/qpcR_datasets/.

using JuliaC

const ROOT = dirname(@__DIR__)
const BUNDLE = joinpath(@__DIR__, "qpcr_fit")

rm(BUNDLE; force = true, recursive = true)
JuliaC.main(
    [
        "--output-exe", "qpcr_fit",
        "--project", ROOT,
        "--bundle", BUNDLE,
        "--trim=safe", "--experimental",
        joinpath(ROOT, "bin", "qpcr_fit.jl"),
    ]
)
cp(joinpath(ROOT, "qpcR_datasets"), joinpath(BUNDLE, "share", "qpcR_datasets"))
cp(joinpath(ROOT, "LICENSE"), joinpath(BUNDLE, "share", "LICENSE"))
println("Executable: ", joinpath(BUNDLE, "bin", Sys.iswindows() ? "qpcr_fit.exe" : "qpcr_fit"))
