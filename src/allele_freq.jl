# allele_freq.jl
# ==============
#
# Compute allele frequencies with BioJulia data types.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

"""
    gene_frequencies(seqcounts::Composition{T}) where T <: Sequence

Compute gene_frequencies from a `BioSequences.Composition` variable that contains
unique sequence counts.
"""
function gene_frequencies(seqcounts::Composition{T}) where T <: Sequence
    n = sum(values(seqcounts))
    frequencies = Dict{T, Float64}()
    @inbounds for (seq, count) in seqcounts
        frequencies[seq] = count / n
    end
    return frequencies
end

"""
    gene_frequencies(iterable)

Compute the gene frequencies for any iterable with an `eltype` which is a
concrete subtype of the abstract `Sequence` type.
"""
function gene_frequencies(iterable)
    return _gene_frequencies(iterable, eltype(iterable), Base.iteratorsize(iterable))
end

# Default for most iterables, throws an error.
_gene_frequencies(iterable, eltype, is) = error("Iterable not supported.")

# Action to take for a sequence iterable which has no known size.
function _gene_frequencies(iterable, ::Type{<:Sequence}, is::Base.SizeUnknown)
    composition = BioSequences.composition(iterable)
    return gene_frequencies(composition)
end

# Action to take for a sequence iterable which has a known size.
# This version computes frequencies directly as n is known in advance.
# The other method first computes composition, and computes frequencies from
# that.
function _gene_frequencies(sequences, ::Type{T}, is::Union{Base.HasLength,Base.HasShape}) where T<:Sequence
    inc = 1 / length(sequences)
    frequencies = Dict{T, Float64}()
    @inbounds for seq in sequences
        old = get(frequencies, seq, 0.0)
        frequencies[seq] = old + inc
    end
    return frequencies
end
