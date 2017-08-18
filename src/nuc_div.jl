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
    nuc_div(m::M, f::V) where {M<:AbstractMatrix{Float64},V<:AbstractVector{Float64}}

Compute nucleotide diversity using a matrix of the number of mutations
between sequence pairs, and a vector of the frequencies of each sequence
in the population.
"""
function nuc_div(m::AbstractMatrix{Float64}, f::AbstractVector{Float64})
    π = 0.0
    @inbounds for i = 1:endof(f), j = (i + 1):endof(f)
        π += m[i, j] * f[i] * f[j]
    end
    return 2 * π
end

"""
    nuc_div(sequences)

Compute nucleotide diversity from any iterable that yields biosequence types.
"""
function nuc_div(sequences)
    frequencies = allele_frequencies(sequences)
    unique_sequences = collect(keys(frequencies))
    n = length(unique_sequences)
    if n < 2
        return 0.0
    else
        mutations = pdist(BioSequences.count_pairwise(Mutated, unique_sequences...))
        return nuc_div(mutations, collect(values(frequencies)))
    end
end
