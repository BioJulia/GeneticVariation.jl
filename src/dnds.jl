# dnds.jl
# =======
#
# dNdS computation.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

const CDN = Union{BioSequences.DNACodon, BioSequences.RNACodon}
const DEFAULT_TRANS = BioSequences.ncbi_trans_table[1]
const CDN_POS_MASKS = (0xFFFFFFFFFFFFFFCF, 0xFFFFFFFFFFFFFFF3, 0xFFFFFFFFFFFFFFFC)
const SITE_PERMUTATIONS = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
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

function pairwise_do!(f::Function, x::Vector{B}, dest::Matrix, opt...) where B <: BioSequence
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

"""
    pairwise_do
"""
function pairwise_do(f::Function, x::Vector{B}, dest::Matrix, opt...) where B <: BioSequence
    n = length(x)
    @assert size(dest) == (n, n) "The size of the dest matrix is not appropriate."
    if n >= 2
        results = Matrix{Tuple{Float64, Float64}}(n, n)
        for i in 1:n
            results[i,i] = 0.0, 0.0
            for j in (i + 1):n
                results[i,j] = results[j,i] = NG86(x[i], x[j], k, code)
            end
        end
        return results
    else
        error("At least two sequences are required.")
    end
end

include("NG86.jl")
