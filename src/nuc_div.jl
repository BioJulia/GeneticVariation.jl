# nuc_div.jl
# ==========
#
# Nucleotide diversity estimation using Nei and Li's 1979 method for sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md


function pdist(d::Matrix{Tuple{Int,Int}})
    o = similar(d, Float64)
    @inbounds for i in eachindex(d)
        o[i] = d[i][1] / d[i][2]
    end
    return o
end

@inline function extract_dict!(dict::Dict{K,V}, keys::Vector{K}, values::Vector{V}) where {K,V}
    @inbounds for (i, pair) in enumerate(dict)
        keys[i] = pair[1]
        values[i] = pair[2]
    end
end

@inline function extract_dict(dict::Dict{K,V}) where {K,V}
    keys = Vector{K}(length(dict))
    values = Vector{V}(length(dict))
    extract_dict!(dict, keys, values)
    return keys, values
end

"""
    NL79(m::M, f::V) where {M<:AbstractMatrix{Float64},V<:AbstractVector{Float64}}

Compute nucleotide diversity using a matrix of the number of mutations
between sequence pairs, and a vector of the frequencies of each sequence
in the population.
"""
function NL79(m::AbstractMatrix{Float64}, f::AbstractVector{Float64})
    π = 0.0
    @inbounds for i = 1:endof(f), j = (i + 1):endof(f)
        π += m[i, j] * f[i] * f[j]
    end
    return 2 * π
end

"""
    NL79(sequences)

Compute nucleotide diversity, as first described by Nei and Li in 1979.

This measure is defined as the average number of nucleotide differences per site
between two DNA sequences in all possible pairs in the sample population, and is
often denoted by the greek letter pi.

`Sequences` should be any iterable that yields biosequence types.

# Examples

```jldoctest
julia> testSeqs = [dna"AAAACTTTTACCCCCGGGGG",
                   dna"AAAACTTTTACCCCCGGGGG",
                   dna"AAAACTTTTACCCCCGGGGG",
                   dna"AAAACTTTTACCCCCGGGGG",
                   dna"AAAAATTTTACCCCCGTGGG",
                   dna"AAAAATTTTACCCCCGTGGG",
                   dna"AAAACTTTTTCCCCCGTAGG",
                   dna"AAAACTTTTTCCCCCGTAGG",
                   dna"AAAAATTTTTCCCCCGGAGG",
                   dna"AAAAATTTTTCCCCCGGAGG"]
10-element Array{BioSequences.BioSequence{BioSequences.DNAAlphabet{4}},1}:
 AAAACTTTTACCCCCGGGGG
 AAAACTTTTACCCCCGGGGG
 AAAACTTTTACCCCCGGGGG
 AAAACTTTTACCCCCGGGGG
 AAAAATTTTACCCCCGTGGG
 AAAAATTTTACCCCCGTGGG
 AAAACTTTTTCCCCCGTAGG
 AAAACTTTTTCCCCCGTAGG
 AAAAATTTTTCCCCCGGAGG
 AAAAATTTTTCCCCCGGAGG

 julia> NL79(testSeqs)
 0.096

```
"""
function NL79(sequences)
    frequencies = gene_frequencies(sequences)
    unique_sequences = collect(keys(frequencies))
    n = length(unique_sequences)
    if n < 2
        return 0.0
    else
        mutations = pdist(BioSequences.count_pairwise(Mutated, unique_sequences...))
        return NL79(mutations, collect(values(frequencies)))
    end
end
