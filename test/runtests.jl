module TestGeneticVariation

using Base.Test

import BioCore.Testing:
    get_bio_fmt_specimens,
    random_seq,
    random_interval
using BioSequences, GeneticVariation
import BufferedStreams: BufferedInputStream
import IntervalTrees: IntervalValue
import YAML

function random_seq(::Type{A}, n::Integer) where A <: Alphabet
    nts = alphabet(A)
    probs = Vector{Float64}(length(nts))
    fill!(probs, 1 / length(nts))
    return BioSequence{A}(random_seq(n, nts, probs))
end

fmtdir = get_bio_fmt_specimens()

include("vcf.jl")
include("bcf.jl")
include("site_counting.jl")
include("minhash.jl")
include("ng86.jl")
include("allele_freq.jl")
include("nuc_div.jl")

end # Module TestGeneticVariation
