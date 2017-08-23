# ng86.jl
# =======
#
# dNdS computation using the NG86 method.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

const CDN = Union{BioSequences.DNACodon, BioSequences.RNACodon}
const DEFAULT_TRANS = BioSequences.ncbi_trans_table[1]

"""
    NG86(x, y, k::Float64 = 1.0, code::GeneticCode)

Compute dNdS statistics, using the Nei and Goborjei 1986 method.

This function requires two iterables `x` and `y`, which yield `DNACodon` or
`RNACodon` type variables. These two types are defined in the BioSequences
package.
"""
function NG86(x, y, k::Float64 = 1.0, code::GeneticCode = DEFAULT_TRANS)
    _NG86(x, y, k, code, eltype(x), eltype(y))
end

"""
    NG86(x::BioSequence{A}, y::BioSequence{A}, k::Float64, code::GeneticCode) where {A <: NucAlphs}

Compute dNdS statistics, using the Nei and Goborjei 1986 method.

This method adds conveinience when working with DNA or RNA sequences, by taking
two sequences, and creating two vectors of aligned codons from them. These two
iterables are then passed into the generic NG86 method.
"""
function NG86(x::BioSequence{A}, y::BioSequence{A}, k::Float64, code::GeneticCode) where {A <: NucAlphs}
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

function _NG86(x, y, k::Float64, code::GeneticCode, xtype::C, ytype::C) where C <: CDN
    # Compute S and N: The expected number of synonymous and nonsynonymous sites.
    S_x, N_x = expected_NG86(NG86, x, k, code)
    S_y, N_y = expected_NG86(NG86, y, k, code)
    S = (S_x + S_y) / 2.0
    N = (N_x + N_y) / 2.0
    # Compute S_d and N_d: The observed number of synonymous and nonsynonymous mutations.
    S_d = N_d = 0.0
    @inbounds for (i, j) in zip(x, y)
        S_i, N_i = observed_NG86(i, j, code)
        S_d += S_i
        N_d += N_i
    end
    # P_s and P_n: The proportion of and synonymous and nonsynonymous differences
    P_s = S_d / S
    P_n = N_d / N
    dN = d_(P_n)
    dS = d_(P_s)
    return dN, dS
end

function _NG86_2(x, y, k::Float64, code::GeneticCode, addone::Bool, xtype::C, ytype::C) where C <: CDN
    S = N = 0.0 # Expected no. of syn and nonsyn sites.
    DS = DN = 0.0 # Observed no. of syn and nonsyn mutations.
    # Iterate over every pair of codons.
    @inbounds for (i, j) in zip(x, y)
        si, ni = expected_NG86(i, k, code)
        sj, nj = expected_NG86(j, k, code)
        S += (si + sj)
        N += (ni + nj)
        DSi, DNi = observed_NG86(i, j, code)
        DS += DSi
        DN += DNi
    end
    S = S / 2.0
    N = N / 2.0
    dN = d_(DN / N)
    dS = addone ? d_((DS + 1) / S) : d_(DS / S)
    return dN, dS
end

"""
    expected_NG86(codon::C, k::Float64, code::GeneticCode)

Enumerate the number of expected synonymous and non-synonymous sites present at
a codon.

Each site may be both partially synonymous and non-synonymous.
"""
function expected_NG86(codon::CDN, k::Float64, code::GeneticCode)
    tsn, tvn = classify_neighbors(codon)
    aa = code[codon]
    S = N = 0.0
    for neighbor in tsn
        if code[neighbor] == AA_Term
            N += 1.0
        elseif code[neighbor] == aa
            S += 1.0
        else
            N += 1.0
        end
    end
    for neighbor in tvn
        if code[neighbor] == AA_Term
            N += k
        elseif code[neighbor] == aa
            S += k
        else
            N += k
        end
    end
    normalization = (N + S) / 3
    return (S / normalization), (N / normalization)
end

function expected_NG86(codons, k::Float64 = 1.0, code::GeneticCode = DEFAULT_TRANS)
    return _expected_NG86(codons, k, code, eltype(codons))
end

function _expected_NG86(codons, k::Float64, code::GeneticCode, etype)
    return error("Iterable not supported.")
end

function _expected_NG86(codons, k::Float64, code::GeneticCode, etype::CDN)
    S = N = 0.0
    @inbounds for codon in codons
        S_i, N_i = _expected_NG86(codon, k, code, etype)
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

"""
    find_differences{C<:CDN}(x::C, y::C)

Identify which sites in two codons are different.
"""
function find_differences(x::C, y::C) where C <: CDN
    positions = Int[1,2,3]
    filter!(i -> x[i] != y[i], positions)
    return positions, length(positions)
end

function observed_NG86(x::C, y::C, code::GeneticCode) where C <: CDN
    if x == y # Early escape, codons are the same, no syn or nonsyn mutations.
        return 0.0, 0.0
    else
        S = N = 0.0

        diff_positions, n_diffs = find_differences(x, y) # Which positions are different.

        if n_diffs == 1

            # One site in the two codons is different. It is obvious and simple
            # then to count whether it is a synonymous or nonsynonymous mutation.
            S, N = classify_mutation(x, y, code)
            return S, N

        elseif n_diffs == 2
            S = N = 0.0

            # For two changes, the number of synonymous and non-synonymous
            # differences per codon, sum to 2, there are two pathways,
            # each possible pathway having two steps.

            # For example, comparing CTA and GTT, the possible pathways are:
            # 1: CTA (L) -> GTA (V) -> GTT (V) : 1 nonsynonymous change and 1 synonymous change.
            # 2: CTA (L) -> CTT (L) -> GTT (V) : 1 nonsynonymous change and 1 synonymous change.

            @inbounds for pos in diff_positions
                bases = collect(x)
                bases[pos] = y[pos]
                temp_cdn = C(bases...)
                # Step 1 of pathway.
                S_i, N_i = classify_mutation(x, temp_cdn, code, 0.5)
                S += S_i
                N += N_i
                # Step 2 of pathway.
                S_i, N_i = classify_mutation(temp_cdn, y, code, 0.5)
                S += S_i
                N += N_i
            end

        elseif n_diffs == 3
            S = N = 0.0

            # For two changes, there are 6 pathways, each with three steps.

            # For example, comparing CTA and GAT, the possible pathways are:
            # 1: CTA (L) -> GTA (V) -> GAA (E) -> GAT (D) : 3 nonsynonymous changes.
            # 2: CTA (L) -> GTA (V) -> GTT (V) -> GAT (D) : 2 nonsynonymous and 1 synonymous change.
            # 3: CTA (L) -> CAA (Q) -> GAA (E) -> GAT (D) : 3 nonsynonymous changes.
            # 4: CTA (L) -> CAA (Q) -> CAT (H) -> GAT (D) : 3 nonsynonymous changes.
            # 5: CTA (L) -> CTT (L) -> GTT (V) -> GAT (D) : 2 nonsynonymous changes and 1 synonymous change.
            # 6: CTA (L) -> CTT (L) -> CAT (H) -> GAT (D) : 2 nonsynonymous changes and 1 synonymous change.

            @inbounds for path in permutations([1,2,3])
                bases = collect(x)
                bases[path[1]] = y[path[1]]
                tmp_cdn_a = C(bases...)
                bases = collect(tmp_cdn_a)
                bases[path[2]] = y[path[2]]
                tmp_cdn_b = C(bases...)
                S_i, N_i = classify_mutation(x, tmp_cdn_a, code, 0.5 / 3)
                S += S_i
                N += N_i
                S_i, N_i = classify_mutation(tmp_cdn_a, tmp_cdn_b, code, 0.5 / 3)
                S += S_i
                N += N_i
                S_i, N_i = classify_mutation(tmp_cdn_b, y, code, 0.5 / 3)
                S += S_i
                N += N_i
            end
        end

        return S, N
    end
end

@inline function d_(p::Float64)
    return - 3 / 4 * log(1 - 4.0 / 3 * p)
end
