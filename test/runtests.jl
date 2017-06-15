module TestGeneticVariation

using Base.Test

using BioSequences, GeneticVariation
import BufferedStreams: BufferedInputStream
import YAML

function get_bio_fmt_specimens()
    path = joinpath(dirname(@__FILE__), "BioFmtSpecimens")
    if !isdir(path)
        run(`git clone --depth 1 https://github.com/BioJulia/BioFmtSpecimens.git $(path)`)
    end
end

include("vcf.jl")


end # Module TestGeneticVariation
