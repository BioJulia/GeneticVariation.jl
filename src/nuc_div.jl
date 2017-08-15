# nuc_div.jl
# ==========
#
# Nucleotide diversity estimation using Nei and Li's 1979 method for sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

@inline function _nalleles(alleles, counts, itersize::T) where {T}
    return sum(values(counts))
end

@inline function _nalleles(alleles, counts, itersize::T) where {T<:Union{Base.HasLength,Base.HasShape}}
    return length(alleles)
end

function allele_frequencies(alleles)
    allele_counts = composition(alleles)
    # Do a bit of dispatch here based on some iterator properties,
    # as not all iterables have a length or size defined. In such a case,
    # we fall back to a less efficient but equally value method.
    n = _nalleles(alleles, allele_counts, Base.iteratorsize(typeof(alleles)))
    frequencies = Dict{eltype(alleles), Float64}()
    @inbounds for (allele, count) in allele_counts
        frequencies[allele] = count / n
    end
    return frequencies
end

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
function nuc_div(m::M, f::V) where {M<:AbstractMatrix{Float64},V<:AbstractVector{Float64}}
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
    mutations = pdist(count_pairwise(Mutated, unique_sequences...))
    return nuc_div(mutations, collect(values(frequencies)))
end
