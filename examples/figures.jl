# Generate the README / documentation figures in docs/src/assets/.
# Usage (from the repository root):
#   julia --project=examples -e 'using Pkg; Pkg.instantiate()'
#   julia --project=examples examples/figures.jl

using CairoMakie
using qPCRModels

const ASSETS = joinpath(@__DIR__, "..", "docs", "src", "assets")
const MODELS = ("Birch" => BirchModel, "Logistic" => LogisticModel, "Richards" => RichardsModel, "MAK2" => MAK2Model)
# Okabe–Ito colors; line styles and markers double the encoding for colour-blind readers
const COLORS = Dict("Birch" => "#0072B2", "Logistic" => "#E69F00", "Richards" => "#009E73", "MAK2" => "#CC79A7")
const LINESTYLES = Dict("Birch" => :solid, "Logistic" => :dash, "Richards" => :dot, "MAK2" => :dashdot)
const MARKERS = Dict("Birch" => :circle, "Logistic" => :rect, "Richards" => :diamond, "MAK2" => :utriangle)

d = qpcRDataset("dil4reps94")
data = d[:, 2:end]
names = d.colnames[2:end]
# column FX_Y: X is the dilution group; F15 amplifies first, so log₁₀(1/X) is the relative template amount
logdil = [-log10(parse(Int, split(split(nm, "_")[1], "F")[2])) for nm in names]
keep = map(!isnoise, eachcol(data))
all(keep) || @warn "$(count(!, keep)) curves show no amplification and are left out"
data, names, logdil = data[:, keep], names[keep], logdil[keep]

fits = Dict(name => [M(Vector(y)) for y in eachcol(data)] for (name, M) in MODELS)
foreach(ms -> foreach(fit!, ms), values(fits))

### Figure 1: the four models fitted to one typical curve

# typical curve: median Birch reduced χ² within the F150 group
group = findall(nm -> startswith(nm, "F150_"), names)
χ²ν = map(reducedchisq, fits["Birch"][group])
j = group[sortperm(χ²ν)[cld(length(group), 2)]]
y = data[:, j]
cycles = 1:length(y)

fig = Figure(; size = (1000, 750), fontsize = 16)
ax1 = Axis(fig[1, 1]; ylabel = "Fluorescence", title = "Curve $(names[j]) (dataset dil4reps94)")
ax2 = Axis(fig[2, 1]; xlabel = "Cycle", ylabel = "Residual")
linkxaxes!(ax1, ax2)
hidexdecorations!(ax1; grid = false, ticks = false)
scatter!(ax1, cycles, y; color = :black, markersize = 9, label = "Data")
hlines!(ax2, 0.0; color = :gray60)
for (name, _) in MODELS
    m = fits[name][j]
    ŷ = get_predicted(m; scaled = false)
    lines!(ax1, 1:length(ŷ), ŷ; color = COLORS[name], linestyle = LINESTYLES[name], linewidth = 2, label = name)
    r = residuals(m; scaled = false)
    scatterlines!(ax2, 1:length(r), r; color = COLORS[name], linestyle = LINESTYLES[name], marker = MARKERS[name], markersize = 8, linewidth = 1.5)
end
nmak2 = fits["MAK2"][j].n
vlines!(ax1, nmak2 + 0.5; color = :gray60, linestyle = :dot)
text!(ax1, nmak2 + 0.5, minimum(y); text = " end of MAK2 window", align = (:left, :bottom), color = :gray40, fontsize = 13)
axislegend(ax1; position = :lt)
rowsize!(fig.layout, 2, Relative(1 / 3))
rowgap!(fig.layout, 8)
save(joinpath(ASSETS, "fit_example.png"), fig; px_per_unit = 2)

### Figure 3: extrapolation of the same curve back to cycle 0, on a log scale

mB, mR = fits["Birch"][j], fits["Richards"][j]
c = 0:length(y)
log10xB = [get_logx0(mB); log.(mB.xmat[:, 2])] ./ log(10)
a, tₘ, ν = get_params(mR)
log10xR = [qPCRModels.logRichardsLogistic(a * (ci - tₘ), ν) for ci in c] ./ log(10)
F_bg, F_max = mB.σ .* mB.linearparams
xdata = (y .- F_bg) ./ F_max # data on the Birch scale; baseline noise is only shown where positive
fig = Figure(; size = (1000, 550), fontsize = 16)
ax3 = Axis(fig[1, 1]; xlabel = "Cycle", ylabel = L"\log_{10}(x)", title = "Curve $(names[j]): extrapolation to cycle 0")
pos = xdata .> 0
scatter!(ax3, cycles[pos], log10.(xdata[pos]); color = :black, markersize = 9, label = "Data (Birch baseline and scale)")
for (name, l) in (("Birch", log10xB), ("Richards", log10xR))
    m = fits[name][j]
    lines!(ax3, c, l; color = COLORS[name], linestyle = LINESTYLES[name], linewidth = 2,
           label = "$name, log₁₀ x₀ = $(round(get_logx0(m) / log(10); digits = 1))")
end
ylims!(ax3, -12, 0.5)
axislegend(ax3; position = :rb)
save(joinpath(ASSETS, "extrapolation_example.png"), fig; px_per_unit = 2)

### Figure 2: log₁₀ x₀ against the dilution, for each model

fig = Figure(; size = (1200, 700), fontsize = 16)
for (k, (name, _)) in enumerate(MODELS)
    ms = fits[name]
    logx0 = map(m -> get_logx0(m) / log(10), ms)
    ok = map(hasconverged, ms)
    cal = calibration(exp10.(logdil), ms) # robust (bisquare) fit of the converged models
    ax = Axis(fig[fldmod1(k, 2)...]; xticks = -5:1:-1, title = name,
              xlabel = L"\log_{10}(1/\mathrm{dilution})", ylabel = L"\log_{10}(x_0)")
    scatter!(ax, logdil[ok], logx0[ok]; color = (COLORS[name], 0.6), marker = MARKERS[name], markersize = 9)
    any(.!ok) && scatter!(ax, logdil[.!ok], logx0[.!ok]; color = :black, marker = :xcross, markersize = 10, label = "Not converged ($(count(!, ok)))")
    ablines!(ax, cal.intercept, cal.slope; color = :black, linestyle = :dash, label = "Slope $(round(cal.slope; digits = 3))\nR² $(round(cal.r2; digits = 3))")
    axislegend(ax; position = (cal.slope > 0 ? :lt : :lb))
end
Label(fig[0, :], "Dataset: dil4reps94 (expected slope: 1)"; fontsize = 20)
colgap!(fig.layout, 16)
rowgap!(fig.layout, 10)
save(joinpath(ASSETS, "x0_dil4reps94.png"), fig; px_per_unit = 2)
