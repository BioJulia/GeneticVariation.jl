# allele_freq.jl
# ==============
#
# Compute allele frequencies with BioJulia data types.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md


"""
    gene_frequencies(iterable)

Compute the gene frequencies for any iterable with an `eltype` which is a
concrete subtype of the abstract `Sequence` type.
"""
function gene_frequencies(iterable)
    return _gene_frequencies(iterable, eltype(iterable))
end

_gene_frequencies(iterable, eltype) = error("Iterable not supported.")

function _gene_frequencies(sequences, ::Type{T}) where T <: Sequence
    seq_counts = BioSequences.composition(sequences)
    n = _ngenes(sequences, seq_counts, Base.iteratorsize(typeof(sequences)))
    frequencies = Dict{T, Float64}()
    @inbounds for (seq, count) in seq_counts
        frequencies[seq] = count / n
    end
    return frequencies
end

@inline function _ngenes(alleles, counts, itersize::T) where {T}
    return sum(values(counts))
end

@inline function _ngenes(alleles, counts, itersize::T) where {T<:Union{Base.HasLength,Base.HasShape}}
    return length(alleles)
end
