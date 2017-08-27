# NG86.jl
# =======
#
# dNdS computation using the NG86 method.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

"""
    dNdS_NG86(x, y, k::Float64 = 1.0, code::GeneticCode)

Compute dNdS statistics, using the Nei and Goborjei 1986 method.

This function requires two iterables `x` and `y`, which yield `DNACodon` or
`RNACodon` type variables. These two types are defined in the BioSequences
package.
"""
function dNdS_NG86(x, y, k::Float64 = 1.0, code::GeneticCode = DEFAULT_TRANS, addone::Bool = false)
    _NG86(x, y, k, code, addone, eltype(x), eltype(y))
end

"""
    dNdS_NG86(x::BioSequence{A}, y::BioSequence{A}, k::Float64, code::GeneticCode) where {A <: NucAlphs}

Compute dNdS statistics, using the Nei and Goborjei 1986 method.

This method adds conveinience when working with DNA or RNA sequences, by taking
two sequences, and creating two vectors of aligned codons from them. These two
iterables are then passed into the generic NG86 method.
"""
function dNdS_NG86(x::BioSequence{A}, y::BioSequence{A}, opt...) where {A <: NucAlphs}
    xcdns, ycdns = aligned_codons(x, y)
    return NG86(xcdns, ycdns, opt...)
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

function _dNdS_NG86(x, y, k::Float64, code::GeneticCode, addone::Bool, xtype::Type{C}, ytype::Type{C}) where C <: CDN
    # Compute S and N: The expected number of synonymous and nonsynonymous sites.
    S_x, N_x = expected_NG86(x, k, code)
    S_y, N_y = expected_NG86(y, k, code)
    S = (S_x + S_y) / 2.0
    N = (N_x + N_y) / 2.0
    # Compute DS and DN: The observed number of synonymous and nonsynonymous mutations.
    DS = ifelse(addone, 1.0, 0.0)
    DN = 0.0
    @inbounds for (i, j) in zip(x, y)
        DS_i, DN_i = observed_NG86(i, j, code)
        DS += DS_i
        DN += DN_i
    end
    # P_s and P_n: The proportion of and synonymous and nonsynonymous differences
    dN = d_(DN / N)
    dS = d_(DS / S)
    return dN, dS
end

function _dNdS_NG86_2(x, y, k::Float64, code::GeneticCode, addone::Bool, xtype::Type{C}, ytype::Type{C}) where C <: CDN
    # Expected no. of syn and nonsyn sites.
    S = N = 0.0
    # Observed no. of syn and nonsyn mutations.
    DS = ifelse(addone, 1.0, 0.0)
    DN = 0.0
    # Iterate over every pair of codons.
    @inbounds for (i, j) in zip(x, y)
        si, ni = S_N_NG86(i, k, code)
        sj, nj = S_N_NG86(j, k, code)
        S += (si + sj)
        N += (ni + nj)
        DSi, DNi = DS_DN_NG86(i, j, code)
        DS += DSi
        DN += DNi
    end
    S = S / 2.0
    N = N / 2.0
    dN = d_(DN / N)
    dS = d_(DS / S)
    return dN, dS
end

"""
    S_N_NG86(codon::C, k::Float64, code::GeneticCode) where {C <: CDN}

Enumerate the number of synonymous (S) and non-synonymous (N) sites in a codon,
using the method used by the Nei and Goborjei (1986).

Returns a tuple where S is the first element and N is the second (S, N). 

Each site in a codon may be both partially synonymous and non-synonymous.
"""
function S_N_NG86(codon::C, k::Float64, code::GeneticCode) where {C <: CDN}
    cdn_bits = UInt64(codon)
    aa = code[codon]
    S = N = 0.0
    for (pos, msk) in enumerate(CDN_POS_MASKS)
        bidx = bitindex(codon, pos)
        @inbounds for base in 0:3
            # Create the neighbor codon.
            neighbor = C((cdn_bits & msk) | (base << bidx))
            if codon == neighbor # Codon created is not a neighbor: should happen 3 times.
                continue
            end
            # See if the mutation is transition or transversion.
            cdn_purine = ispurine(codon[pos])
            neighbor_purine = ispurine(neighbor[pos])
            istransition = (cdn_purine && neighbor_purine) || (!cdn_purine && !neighbor_purine)
            # See if the protein changes between codon and neighbor, and update
            # N and S counts accordingly.
            inc = ifelse(istransition, 1.0, k)
            neighbor_aa = code[neighbor]
            if neighbor_aa == AA_Term
                N += inc
            elseif neighbor_aa == aa
                S += inc
            else
                N += inc
            end
        end
    end
    normalization = (N + S) / 3
    return (S / normalization), (N / normalization)
end

@inline bitindex(x::Kmer{T,K}, i::Integer) where {T,K} = 2 * (K - i)

function S_N_NG86(codons, k::Float64 = 1.0, code::GeneticCode = DEFAULT_TRANS)
    return _expected_NG86(codons, k, code, eltype(codons))
end

function _S_N_NG86(codons, k::Float64, code::GeneticCode, etype)
    return error("Iterable not supported.")
end

function _S_N_NG86(codons, k::Float64, code::GeneticCode, etype::Type{C}) where C <: CDN
    S = N = 0.0
    @inbounds for codon in codons
        S_i, N_i = expected_NG86(codon, k, code)
        S += S_i
        N += N_i
    end
    return S, N
end

@inline function classify_mutation(x::C, y::C, code::GeneticCode, weight::Float64 = 1.0) where C <: CDN
    if code[x] == code[y]
        # Synonymous mutation.
        return weight, 0.0
    else
        # non-synonymous mutation.
        return 0.0, weight
    end
end

function find_differences(x::C, y::C) where C <: CDN
    diffs = 0x00
    @inbounds for pos in 1:3
        diffs = (x[pos] != y[pos]) | (diffs << 1)
    end
    return diffs, count_ones(diffs)
end

function splice_into(x::C, y::C, pos::Integer) where {C <: CDN}
    mask = UInt64(3) << bitindex(x, pos)
    return C((UInt64(x) & ~mask) | (UInt64(y) & mask))
end

"""
    DS_DN_NG86(x::C, y::C, code::GeneticCode) where C <: CDN

Compute the number of synonymous (DS) and non-synonymous (DN) mutations between
two codons, using the all paths method used by the Nei and Goborjei (1986).
"""
function DS_DN_NG86(x::C, y::C, code::GeneticCode) where C <: CDN
    if x == y # Early escape, codons are the same, no syn or nonsyn mutations.
        return 0.0, 0.0
    else
        diff_positions, n_diffs = find_differences(x, y) # Which positions are different.
        if n_diffs == 1
            # One site in the two codons is different. It is obvious and simple
            # then to count whether it is a synonymous or nonsynonymous mutation.
            DS, DN = classify_mutation(x, y, code)
            return DS, DN
        elseif n_diffs == 2
            DS = DN = 0.0
            # For two changes, the number of synonymous and non-synonymous
            # differences per codon, sum to 2, there are two pathways,
            # each possible pathway having two steps.
            # For example, comparing CTA and GTT, the possible pathways are:
            # 1: CTA (L) -> GTA (V) -> GTT (V) : 1 nonsynonymous change and 1 synonymous change.
            # 2: CTA (L) -> CTT (L) -> GTT (V) : 1 nonsynonymous change and 1 synonymous change.
            @inbounds for pos in 1:3
                if (diff_positions & (0x01 << (pos - 1))) > 0x00
                    temp_cdn = splice_into(x, y, pos)
                    # Step 1 of pathway.
                    DS_i, DN_i = classify_mutation(x, temp_cdn, code, 0.5)
                    DS += DS_i
                    DN += DN_i
                    # Step 2 of pathway.
                    DS_i, DN_i = classify_mutation(temp_cdn, y, code, 0.5)
                    DS += DS_i
                    DN += DN_i
                end
            end
        elseif n_diffs == 3
            DS = DN = 0.0
            # For two changes, there are 6 pathways, each with three steps.
            # For example, comparing CTA and GAT, the possible pathways are:
            # 1: CTA (L) -> GTA (V) -> GAA (E) -> GAT (D) : 3 nonsynonymous changes.
            # 2: CTA (L) -> GTA (V) -> GTT (V) -> GAT (D) : 2 nonsynonymous and 1 synonymous change.
            # 3: CTA (L) -> CAA (Q) -> GAA (E) -> GAT (D) : 3 nonsynonymous changes.
            # 4: CTA (L) -> CAA (Q) -> CAT (H) -> GAT (D) : 3 nonsynonymous changes.
            # 5: CTA (L) -> CTT (L) -> GTT (V) -> GAT (D) : 2 nonsynonymous changes and 1 synonymous change.
            # 6: CTA (L) -> CTT (L) -> CAT (H) -> GAT (D) : 2 nonsynonymous changes and 1 synonymous change.
            @inbounds for path in SITE_PERMUTATIONS
                tmp_cdn_a = splice_into(x, y, path[1])
                tmp_cdn_b = splice_into(tmp_cdn_a, y, path[2])
                DS_i, DN_i = classify_mutation(x, tmp_cdn_a, code, 0.5 / 3)
                DS += DS_i
                DN += DN_i
                DS_i, DN_i = classify_mutation(tmp_cdn_a, tmp_cdn_b, code, 0.5 / 3)
                DS += DS_i
                DN += DN_i
                DS_i, DN_i = classify_mutation(tmp_cdn_b, y, code, 0.5 / 3)
                DS += DS_i
                DN += DN_i
            end
        end
        return DS, DN
    end
end

@inline function d_(p::Float64)
    return - 3 / 4 * log(1 - 4.0 / 3 * p)
end