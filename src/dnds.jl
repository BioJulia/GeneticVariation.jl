# dnds.jl
# =======
#
# dNdS computation.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

const CDN = Union{BioSequences.DNACodon, BioSequences.RNACodon}
const DEFAULT_TRANS = BioSequences.ncbi_trans_table[1]
#=
struct AlignedCodons{T<:NucAlphs}
    x::BioSequence{T}
    y::BioSequence{T}
end

@inline start(ac::AlignedCodons)::Int = 1
@inline function next(ac::AlignedCodons{T}, state::Int) where T<: NucAlphs
    cdnx, okx = BioSequences.extract_kmer_impl(x, state, 3)
    cdny, oky = BioSequences.extract_kmer_impl(y, state, 3)
    state += 3
    if okx && oky
        return (cdnx, cdny), state
    else
        return next(ac, state)
    end
end
@inline function done(ac::AlignedCodons, state::Int)::Bool
    return state + 2 > min(endof(ac.x), endof(ac.y))
end
=#
"""
    aligned_codons(x::BioSequence{T}, y::BioSequence{T}, start::Int = 1) where T <: NucAlphs

Create two
"""
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

@inline bases(::Type{BioSequences.DNACodon}) = BioSequences.ACGT
@inline bases(::Type{BioSequences.RNACodon}) = BioSequences.ACGU

"""
    classify_neighbor(codon::DNACodon)

Compute and classify the neighbors of a given `codon` as either a transitio
neighbor, or a transversion neighbor.
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

function pairwise_do!(f::Function, x::Vector{B}, dest::Matrix, opt...) where B<:BioSequence
    n = length(x)
    @assert size(dest) == (n, n) "The size of the dest matrix is not appropriate."
    @assert l >= 2 "Not enough sequences."
    for i in 1:n
        dest[i,i] = zero(eltype(dest))
        for j in (i + 1):n
            dest[i,j] = dest[j,i] = f(x[i], x[j], opt...)
        end
    end
    return dest
end

include("NG86.jl")
