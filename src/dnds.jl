# dnds.jl
# =======
#
# dNdS computation.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

const CDN = Union{BioSequences.DNACodon, BioSequences.RNACodon}
const DEFAULT_TRANS = BioSequences.ncbi_trans_table[1]

@inline bases(::Type{BioSequences.DNACodon}) = BioSequences.ACGT
@inline bases(::Type{BioSequences.RNACodon}) = BioSequences.ACGU

function aligned_codons(x::BioSequence{T}, y::BioSequence{T}, start::Int = 1) where T <: NucAlphs
    xcdns = Vector{Kmer{eltype(T), 3}}()
    ycdns = Vector{Kmer{eltype(T), 3}}()
    pos = start
    while pos + 2 ≤ min(endof(x), endof(y))
        cdnx, okx = BioSequences.extract_kmer_impl(x, pos, 3)
        cdny, oky = BioSequences.extract_kmer_impl(y, pos, 3)
        if okx && oky
            push!(xcdns, convert(Kmer{eltype(T), 3}, cdnx))
            push!(ycdns, convert(Kmer{eltype(T), 3}, cdny))
        end
        pos += 3
    end
    return xcdns, ycdns
end

"""
    classify_neighbor(codon::DNACodon)

Computes and classifies the neighbors of a given `codon` as either a
transition neighbor, or a transversion neighbor.
"""
function classify_neighbors(codon::C) where C <: CDN
    tsn = Vector{C}()
    tvn = Vector{C}()
    codon_bases = collect(codon)
    @inbounds for n in 1:3
        i = codon_bases[n]
        for j in bases(C)
            if i != j
                thiscdn = copy(codon_bases)
                thiscdn[n] = j
                ipur = ispurine(i)
                jpur = ispurine(j)
                topush = ifelse((ipur && jpur) || (!ipur && !jpur), tsn, tvn)
                push!(topush, C(thiscdn...))
            end
        end
    end
    return tsn, tvn
end

function dNdS(::Type{T}, x::BioSequence{A}, y::BioSequence{A}, opt...) where {T, A <: NucAlphs}
    return dNdS(T, aligned_codons(x, y)..., opt...)
end

function pairwise_dNdS(::Type{T}, x::Vector{V}, opt...) where {T, V}
    n = length(x)
    @assert n >= 2 "At least two sequences are required."
    results = Matrix{Tuple{Float64, Float64}}(n, n)
    for i in 1:n
        results[i,i] = 0.0, 0.0
        for j in (i + 1):n
            results[i,j] = results[j,i] = dNdS(T, x[i], x[j], opt...)
        end
    end
    return results
end

include("ng86.jl")
