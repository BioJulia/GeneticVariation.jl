# site_counting.jl
# ================
#
# Extension to the site counting framework for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

@compat abstract type Mutation <: Position end

"""
A `Conserved` site describes a site where two aligned nucleotides are definately
conserved. By definately conserved this means that the symbols of the site are
non-ambiguity symbols, and they are the same symbol.
"""
immutable Conserved <: Mutation end

"""
A `Mutated` site describes a site where two aligned nucleotides are definately
mutated. By definately mutated this means that the symbols of the site are
non-ambiguity symbols, and they are not the same symbol.
"""
immutable Mutated <: Mutation end

"""
A `Transition` site describes a site where two aligned nucleotides are definately
mutated, and the type of mutation is a transition mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transition mutation: i.e. A<->G, or C<->T.
"""
immutable Transition <: Mutation end

"""
A `Transversion` site describes a site where two aligned nucleotides are
definately mutated, and the type of mutation is a transversion mutation.
In other words, the symbols must not be ambiguity symbols, and they must
be different such that they constitute a transversion mutation: i.e. A<->C,
A<->T, G<->T, G<->C.
"""
immutable Transversion <: Mutation end


# Masking functions
# -----------------

@inline function nibble_mask(::Type{Certain}, x::UInt64)
    return nibble_mask(enumerate_nibbles(x), 0x1111111111111111)
end

@inline function nibble_mask(::Type{Certain}, a::UInt64, b::UInt64)
    return nibble_mask(Certain, a) & nibble_mask(Certain, b)
end


# BioSequences.count_sites_bitpar extension
# -----------------------------------------

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin

        # Counter types
        @inline BioSequences.bp_counter_type{M<:Mutation}(::Type{M}, ::Type{$A{4}}) = Tuple{Int, Int}
        @inline BioSequences.bp_start_counter{M<:Mutation}(::Type{M}, ::Type{$A{4}}) = Int(0), Int(0)

        # Conserved
        @inline bp_chunk_count(::Type{Conserved}, ::Type{$A{2}}, a::UInt64, b::UInt64) = bp_chunk_count(Match, $A{2}, a, b)

        @inline function bp_chunk_count(::Type{Conserved}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            k, c = bp_chunk_count(Mutated, $A{4}, a, b)
            return c - k, c
        end

        @inline correct_emptyspace(::Type{Conserved}, ::Type{$A{2}}) = true

        # Mutated
        @inline bp_chunk_count(::Type{Mutated}, ::Type{$A{2}}, a::UInt64, b::UInt64) = bp_chunk_count(Mismatch, $A{2}, a, b)

        @inline function bp_chunk_count(::Type{Mutated}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            m = nibble_mask(Certain, a, b)
            d = (a $ b) & m
            return count_nonzero_nibbles(d), count_one_nibbles(m)
        end

        # TODO Finish transitions and transversions
        # Transition
        @inline function bp_chunk_count(::Type{Transition}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            m = nibble_mask(Certain, a, b)
            d = (a $ b) & m
            remcount = count_zero_nibbles(d) # Count of the identical nucleotides.
            tve = count_one_nibbles(nibble_mask(0x9999999999999999, d)) # Count the 1001 transversion edge case.
            d &= (d >> 1)
            d &= 0x7777777777777777
            tvc = count_ones(d)
            tve += (16 - tvc - remcount)
        end


        # Transversion
        # TODO

    end
end

@inline BioSequences.bp_update_counter(acc::Tuple{Int,Int}, up::Tuple{Int,Int}) = acc[1] + up[1], acc[2] + up[2]
