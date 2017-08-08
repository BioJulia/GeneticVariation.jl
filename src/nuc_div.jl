# nuc_div.jl
# ==========
#
# Nucleotide diversity estimation using Nei and Li's 1979 method.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

"""
    count_unique_sequences(seqs::Vararg{BioSequence{A},N}) where {A <: Alphabet, N}

Count the number of unique sequences in a set of sequences.
"""
@inline function count_unique_sequences(seqs::Vararg{BioSequence{A},N}) where {A <: Alphabet, N}
    counts = Dict{BioSequence{A}, Int}()
    @inbounds for seq in seqs
        current = get(counts, seq, 0)
        counts[seq] = current + 1
    end
    return counts
end

@inline function count_unique_sequences(seqs::Vector{BioSequence{A},N}) where {A <: Alphabet, N}
    count_unique_sequences(seqs...)
end

@inline function allele_frequencies(seqs::Vararg{BioSequence{A},N}) where {A <: NucAlphs, N}
    seq_counts = count_unique_sequences(seqs)
    frequencies = Dict{BioSequence{A}, Float64}()
    @inbounds for seq, count in seq_counts
        frequencies[seq] = count / N
    end
    return frequencies
end

#=
@inline function nuc_div(seqs::Vararg{BioSequence{A},N}) where {A<:NucAlphs, N}
    scounts = count_unique_sequences(seqs)

end
=#


"""
    nuc_div(d::Matrix{Int}, f::Vector{Float64})

Compute nucleotide diversity using a matrix of the number of mutations
between sequence pairs, and a vector of the frequencies of each sequence
in the population.
"""
@inline function nuc_div(d::Matrix{Int}, f::Vector{Float64})
    tot = 0.0
    @inbounds @simd for i = 1:endof(f), j = (i + 1):endof(f)
        tot += d[i, j] * f[i] * f[j]
    end
    return 2 * tot
end

"""
    nuc_div(d::Matrix{Int}, f::Vector{Float64})

Compute nucleotide diversity using only a matrix of the number of mutations
between sequence pairs.

This assumes each sequence was present in equal frequencies
in the population.
"""
@inline function nuc_div(d::Matrix{Int})
    dsize = size(d)[2]
    return nuc_div(d, fill(1 / dsize, dsize))
end

function nuc_div(s::Vector{DNASequence})
    nuc_div(extract_first())
end

@inline function extract_first(d::Matrix{Tuple{Int,Int}})
    o = similar(d, Int)
    @inbounds @simd for i in eachindex(d)
        o[i] = d[i][1]
    end
    return o
end
