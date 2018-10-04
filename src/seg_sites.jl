# seg_sites.jl
# =============
#
# Identify and count segregating sites with BioJulia data types.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

"""
`Segregating` sites are positions which show differences (polymorphisms)
between related genes in a sequence alignment (are not conserved).
Segregating sites include conservative, semi-conservative and non-conservative
mutations.
"""
struct Segregating <: Position end

function Base.count(::Type{Segregating}, seqs::Vector{T}) where T <: Sequence
    cont = true
    i = 0
    count = 0
    @inbounds refseq = seqs[1]
    while cont
        i += 1
        @inbounds refelem = refseq[i]
        @inbounds for j in 2:lastindex(seqs)
            seq = seqs[j]
            cont = ifelse(lastindex(seq) == i, false, true)
            if refelem != seq[i]
                count += 1
                break
            end
        end
    end
    return count, i
end
