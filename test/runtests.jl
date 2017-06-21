module TestGeneticVariation

using Base.Test

using BioSequences, GeneticVariation
import BufferedStreams: BufferedInputStream
import IntervalTrees: IntervalValue
import YAML

function get_bio_fmt_specimens()
    path = joinpath(dirname(@__FILE__), "BioFmtSpecimens")
    if !isdir(path)
        run(`git clone --depth 1 https://github.com/BioJulia/BioFmtSpecimens.git $(path)`)
    end
end

function random_seq(n::Integer, nts, probs)
    cumprobs = cumsum(probs)
    x = Vector{Char}(n)
    for i in 1:n
        x[i] = nts[searchsorted(cumprobs, rand()).start]
    end
    return convert(String, x)
end

function random_seq{A<:Alphabet}(::Type{A}, n::Integer)
    nts = alphabet(A)
    probs = Vector{Float64}(length(nts))
    fill!(probs, 1 / length(nts))
    return BioSequence{A}(random_seq(n, nts, probs))
end

function random_interval(minstart, maxstop)
    start = rand(minstart:maxstop)
    return start:rand(start:maxstop)
end

include("vcf.jl")
include("bcf.jl")
include("site_counting.jl")


end # Module TestGeneticVariation
