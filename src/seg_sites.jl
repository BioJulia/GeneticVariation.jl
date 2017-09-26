# seg_sites.jl
# =============
#
# Identify and count segregating sites with BioJulia data types.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

"""
    nsegregating(seqs::Vector{T}) where T <: Sequence

Compute the number of segregating sites in a vector of sequences.
The set of sequences are assumed to be aligned, starting with their first
position.
"""
function nsegregating(seqs::Vector{T}) where T <: Sequence
    cont = true
    i = 0
    count = 0
    @inbounds refseq = seqs[1]
    while cont
        i += 1
        @inbounds refelem = refseq[i]
        @inbounds for j in 2:endof(seqs)
            seq = seqs[j]
            cont = ifelse(endof(seq) == i, false, true)
            if refelem != seq[i]
                count += 1
                break
            end
        end
    end
    return count, i
end
