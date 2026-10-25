# Dataset loading and qpcR dataset container

"""
    get_dataset_folder()

Return absolute filepath to the built-in `qpcR_datasets/` directory associated with `qPCRModels`.
In a standalone executable built with juliac, the package directory is the one used at build time,
so the copy in the bundle (`share/qpcR_datasets/`, next to `bin/`) is used instead.
"""
function get_dataset_folder()
    pkg_path = pkgdir(@__MODULE__)
    folder = isnothing(pkg_path) ? "qpcR_datasets" : joinpath(pkg_path, "qpcR_datasets")
    bundled = normpath(joinpath(Sys.BINDIR, "..", "share", "qpcR_datasets"))
    return !isdir(folder) && isdir(bundled) ? bundled : folder
end

"""
    DatasetDetails

Dictionary storing metadata annotations for all built-in benchmark qPCR datasets (e.g. `batsch1`, `boggy`, `karlen1`, `testdat`).
Contains metadata on experimental setup, column annotations, hardware/instrumentation, and assay details.
"""
const DatasetDetails = Dict(
    "batsch1" => (Setup = "Five 4-fold dilutions with 3 replicates.",
                Annotation = "FX.Y (X = dilution number, Y = replicate number).", 
                Hardware = "Lightcycler 1.0 (Roche).",
                Details = "Primers for rat SLC6A14, Taqman probes."),
    "batsch2" => (Setup = "Five 4-fold dilutions with 3 replicates.",
                Annotation = "FX.Y (X = dilution number, Y = replicate number).", 
                Hardware = "Lightcycler 1.0 (Roche).",
                Details = "Primers for human SLC22A13, Taqman probes."),
    "batsch3" => (Setup = "Five 4-fold dilutions with 3 replicates.",
                Annotation = "FX.Y (X = dilution number, Y = replicate number).", 
                Hardware = "Lightcycler 1.0 (Roche).",
                Details = "Primers for pig EMT, Taqman probes."),
    "batsch4" => (Setup = "Five 4-fold dilutions with 3 replicates.",
                Annotation = "FX.Y (X = dilution number, Y = replicate number).", 
                Hardware = "Lightcycler 1.0 (Roche).",
                Details = "Primers for chicken ETT, SybrGreen."),
    "batsch5" => (Setup = "Five 4-fold dilutions with 3 replicates.",
                Annotation = "FX.Y (X = dilution number, Y = replicate number).", 
                Hardware = "Lightcycler 1.0 (Roche).",
                Details = "Primers for human GAPDH, SybrGreen."),
    "boggy" => (Setup = "Six 10-fold dilutions with 2 replicates.",
                Annotation = "FX.Y (X = dilution number, Y = replicate number).", 
                Hardware = "Chromo4 (BioRad).",
                Details = "Primers for a synthetic template, consisting of a secondary structure-optimized random sequence (129 bp), Syto-13 dye."),
    "competimer" => (Setup = "7 concentrations of inhibitor, six 4-fold dilutions, 3 replicates.",
                Annotation = "X_Y_Z (X = inhibitor concentration, Y = dilution number, Z = replicate number).\n
                    X: % competimer\n
                    A 0%, B 5%, C 10%, D 20%, E 30%, F 40%, G 50%.\n
                    Y: dilution factor (-fold)\n
                    A 64, B 16, C 4, D 1, E 0.25, F 0.0625, G NTC.", 
                Hardware = "Lightcycler 480 (Roche).",
                Details = "Primers for human AluSx repeats, competitive primers, SybrGreen I dye. NTCs are included."),
    "dil4reps94" => (Setup = "Four 10-fold dilutions with 94 replicates.",
                Annotation = "FX_Y (X = copy number, Y = replicate number).", 
                Hardware = "CFX384 (BioRad).",
                Details = "Primers for the human MYCN gene, synthetic MYCN oligo used as template, SybrGreen I dye. NTCs were removed."),
    "dyemelt" => (Setup = "Melting curves of a 4-plex qPCR with different fluorescence dyes.",
                Annotation = "T0 (Temperature), EvaGreen, T1 (Temperature), SybrGreen.I, T2 (Temperature), Syto13.", 
                Hardware = "Lightcycler 1.0 (Roche).",
                Details = "A melting curve analysis of a 4-plex real-time PCR on genomic DNA with AZF deletion-specific primers. The dyes used were EvaGreen, SybrGreen I and Syto-13."),
    "guescini1" => (Setup = "Seven 10-fold dilutions with 12 replicates.",
                Annotation = "FX.Y (X = dilution number, Y = replicate number).", 
                Hardware = "Lightcycler 480 (Roche).",
                Details = "Primers for NADH dehydrogenase 1, SybrGreen I dye, data is background subtracted."),
    "guescini2" => (Setup = "Five decreasing steps of PCR mix with 12 replicates.",
                Annotation = "FX.Y (X = dilution number, Y = replicate number).", 
                Hardware = "Lightcycler 480 (Roche).",
                Details = "Primers for NADH dehydrogenase 1, SybrGreen I dye, data is background subtracted."),
    "htPCR" => (Setup = "High throughput experiment containing 8858 runs from a 95 x 96 PCR grid.",
                Annotation = "PX.Y (X = plate number, Y = well number).", 
                Hardware = "Biomark HD (Fluidigm).",
                Details = "Proprietary primers, EvaGreen dye, data is ROX normalized."),
    "karlen1" => (Setup = "4 dilutions (1-, 10-, 50-, 100-fold) with 5 replicates in 4 samples.",
                Annotation = "FX.Y.Z (X = sample number, Y = dilution number, Z = replicate number).", 
                Hardware = "ABI Prism 7700 (Applied Biosystems).",
                Details = "Primers for Caveolin, SybrGreen I dye, data is background subtracted."),
    "karlen2" => (Setup = "5 dilutions (1-, 10-, 50-, 100, 1000-fold) with 5 (4 for 1000-fold dilution) replicates in 4 samples.",
                Annotation = "FX.Y.Z (X = sample number, Y = dilution number, Z = replicate number).", 
                Hardware = "ABI Prism 7700 (Applied Biosystems).",
                Details = "Primers for Fibronectin, SybrGreen I dye, data is background subtracted."),
    "karlen3" => (Setup = "5 dilutions (1-, 10-, 50-, 100, 1000-fold) with 5 (4 for 1000-fold dilution) replicates in 4 samples.",
                Annotation = "FX.Y.Z (X = sample number, Y = dilution number, Z = replicate number).", 
                Hardware = "ABI Prism 7700 (Applied Biosystems).",
                Details = "Primers for L27 (karlen3), SybrGreen I dye, data is background subtracted."),
    "lievens1" => (Setup = "Five 5-fold dilutions with 18 replicates.",
                Annotation = "SX.Y (X = dilution number, Y = replicate number).", 
                Hardware = "ABI7300 (ABI) or Biorad IQ5 (Biorad).",
                Details = "Primers for the soybean lectin endogene Le1, SybrGreen I dye."),
    "lievens2" => (Setup = "Five different concentrations of isopropanol (2.5%, 0.5%, 0.1%, 0.02% and 0.004% (v/v)) with 18 replicates.",
                Annotation = "SX.Y (X = concentration step, Y = replicate number).", 
                Hardware = "ABI7300 (ABI) or Biorad IQ5 (Biorad).",
                Details = "Primers for the soybean lectin endogene Le1, SybrGreen I dye."),
    "lievens3" => (Setup = "Five different amounts of tannic acid per reaction (5 ng, 1 ng, 0.2 ng, 0.04 ng and 0.008 ng) and 18 replicates.",
                Annotation = "SX.Y (X = concentration step, Y = replicate number).", 
                Hardware = "ABI7300 (ABI) or Biorad IQ5 (Biorad).",
                Details = "Primers for the soybean lectin endogene Le1, SybrGreen I dye."),
    "reps" => (Setup = "Seven 10-fold dilutions with 4 replicates.",
                Annotation = "FX.Y (X = dilution number, Y = replicate number).", 
                Hardware = "Lightcycler 1.0 (Roche).",
                Details = "Primers for the S27a housekeeping gene, SybrGreen I dye."),
    "reps2" => (Setup = "Five 4-fold dilutions with 3 replicates, 2 different cDNAs.",
                Annotation = "FX.Y.Z (X = cDNA number, Y = dilution number, Z = replicate number).", 
                Hardware = "MXPro3000P (Stratagene).",
                Details = "Primers for the S27a housekeeping gene, SybrGreen I dye."),
    "reps3" => (Setup = "Seven 4-fold dilutions with 3 replicates.",
                Annotation = "FX.Y (X = dilution number, Y = replicate number).", 
                Hardware = "MXPro3000P (Stratagene).",
                Details = "Primers for the S27a housekeeping gene, SybrGreen I dye, ROX-normalized."),
    "reps384" => (Setup = "A data frame with 379 replicate runs of a 384 microtiter plate.",
                Annotation = "A_A_X (X = replicate number).", 
                Hardware = "CFX384 (BioRad).",
                Details = "Primers for the human MYCN gene, synthetic MYCN oligo used as template (15000 copies), SybrGreen I dye. NTCs were removed."),
    "rutledge" => (Setup = "Six 10-fold dilutions with 4 replicates in 5 individual batches.",
                Annotation = "X.RY.Z (X = dilution number, Y = batch number, Z = replicate number).", 
                Hardware = "Opticon 2 (MJ Research).",
                Details = "Primers for a 102 bp amplicon, SybrGreen I dye, data is background subtracted."),
    "testdat" => (Setup = "Six 10-fold dilutions with 4 replicates.",
                Annotation = "FX.Y (X = dilution number, Y = replicate number).", 
                Hardware = "Lightcycler 1.0 (Roche).",
                Details = "Same as reps, but each FX.3 has noisy data which fails to fit with the l5 model, each FX.4 passes fitting but fails in sigmoidal structure detection by KOD. Used for evaluating quality checking methods."),
    "vermeulen1" => (Setup = "A subset of the first 20 samples for each of 64 genes.",
                Annotation = "X.Y (X = gene name, Y = sample name).", 
                Hardware = "Lightcycler 480 (Roche).",
                Details = "Primers for AHCY, AKR1C1, ALUsq(Eurogentec), ARHGEF7, BIRC5, CAMTA1, CAMTA2, CD44, CDCA5, CDH5, CDKN3, CLSTN1, CPSG3, DDC, DPYSL3, ECEL1, ELAVL4, EPB41L3, EPHA5, EPN2, FYN, GNB1, HIVEP2, HMBS, HPRT1, IGSF4, INPP1, MAP2K4, MAP7, MAPT, MCM2, MRPL3, MTSS1, MYCN(4), NHLH2, NM23A, NRCAM, NTRK1, ODC1, PAICS, PDE4DIP, PIK3R1, PLAGL1, PLAT, PMP22, PRAME, PRDM2, PRKACB, PRKCZ, PTN, PTPRF, PTPRH, PTPRN2, QPCT, SCG2, SDHA(1), SLC25A5, SLC6A8, SNAPC1, TNFRSF, TYMS, UBC(2), ULK2 and WSB1. SybrGreen I dye."),
    "vermeulen2" => (Setup = "vermeulen1 corresponding dilution data for all 64 genes with five 10-fold dilutions and 3 replicates.",
                Annotation = "X.STD_Y.Z (X = gene name, Y = copy number, Z = replicate number).", 
                Hardware = "Lightcycler 480 (Roche).",
                Details = "Primers for AHCY, AKR1C1, ALUsq(Eurogentec), ARHGEF7, BIRC5, CAMTA1, CAMTA2, CD44, CDCA5, CDH5, CDKN3, CLSTN1, CPSG3, DDC, DPYSL3, ECEL1, ELAVL4, EPB41L3, EPHA5, EPN2, FYN, GNB1, HIVEP2, HMBS, HPRT1, IGSF4, INPP1, MAP2K4, MAP7, MAPT, MCM2, MRPL3, MTSS1, MYCN(4), NHLH2, NM23A, NRCAM, NTRK1, ODC1, PAICS, PDE4DIP, PIK3R1, PLAGL1, PLAT, PMP22, PRAME, PRDM2, PRKACB, PRKCZ, PTN, PTPRF, PTPRH, PTPRN2, QPCT, SCG2, SDHA(1), SLC25A5, SLC6A8, SNAPC1, TNFRSF, TYMS, UBC(2), ULK2 and WSB1. SybrGreen I dye.")
)

"""
    qpcRDataset{T<:AbstractFloat} <: AbstractMatrix{T}

Container holding a benchmark qPCR dataset loaded from tab-delimited text files.

Subtypes `AbstractMatrix{T}` so data columns (cycles, fluorescence curves) can be indexed directly (`ds[:, 2]`).

# Fields
- `details`: `@NamedTuple{Setup, Annotation, Hardware, Details}` experimental information.
- `colnames::Vector{String}`: Names of data columns (e.g. `["Cycles", "F1.1", "F1.2", ...]`).
- `data::Matrix{T}`: Numeric fluorescence data matrix (rows = cycles, cols = samples/curves).

# Examples
```julia
ds = qpcRDataset("batsch1")
cycle_data = ds[:, 1]
curve1 = ds[:, 2]
```
"""
struct qpcRDataset{T<:AbstractFloat} <: AbstractMatrix{T}
    details::@NamedTuple{Setup::String, Annotation::String, Hardware::String, Details::String}
    colnames::Vector{String}
    data::Matrix{T}

    function qpcRDataset(::T, name::String) where {T<:AbstractFloat}
        filepath = joinpath(get_dataset_folder(), name * ".txt")
        txt = split(read(filepath, String), "\n"; keepempty=false) # readlines(filepath) cannot be trimmed by juliac
        colnames = String.(split(txt[1], "\t"; keepempty=false))
        data = zeros(T, length(txt) - 1, length(colnames))
        for (i, lᵢ) in Iterators.enumerate(Iterators.drop(txt, 1))
            data[i, :] .= parse.(T, split(lᵢ, "\t"; keepempty=false))
        end
        return new{T}(DatasetDetails[name], colnames, data)
    end
end
qpcRDataset(name::String) = qpcRDataset(zero(Float64), name)
Base.size(d::qpcRDataset) = size(d.data)
Base.IndexStyle(::qpcRDataset) = IndexCartesian()
Base.getindex(d::qpcRDataset, I...) = getindex(d.data, I...)

"""
    load_qpcR_dataset([Type{T}=Float64])

Load all available built-in qPCR benchmark datasets into a dictionary mapping dataset names to `qpcRDataset{T}`.

# Arguments
- `T::Type{<:AbstractFloat}`: Floating point type (default: `Float64`).

# Returns
- `Dict{String, qpcRDataset{T}}`: Dictionary of loaded datasets.
"""
function load_qpcR_dataset(::Type{T}) where {T <: AbstractFloat}
    qpcRDatasets = Dict{String, qpcRDataset{T}}()
    for name in keys(DatasetDetails)
        qpcRDatasets[name] = qpcRDataset(zero(T), name)
    end
    return qpcRDatasets
end
load_qpcR_dataset() = load_qpcR_dataset(Float64)
