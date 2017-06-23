# proportion.jl
# =============
#
# Types and methods for computing p-distances.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

"""
Proportion{T} is a distance which is the count of the mutations of type T that
exist between the two biological sequences, divided by the number of valid sites
examined (sites which don't have gap or ambiguous symbols).
In other words this so called p-distance is simply the proportion of sites
between each pair of sequences, that are mutated (again where T determines
what kind of mutation).
"""
immutable Proportion{T} <: UncorrectedDistance end

@inline function distance{T}(::Type{Proportion{T}}, x::Sequence, y::Sequence)
    cs = count(T, x, y)
    return _pcorrection(Proportion{T}, cs, x, y)
end

@inline function pdistance{T}(::Type{T}, x, y)
    return distance(Proportion{T}, x, y)
end

@inline function _pcorrection(k::Int, x::Sequence, y::Sequence)
    return k / min(length(x), length(y))
end

@inline function _pcorrection(cs::Tuple{Int,Int}, x::Sequence, y::Sequence)
    return cs[1] / cs[2]
end
