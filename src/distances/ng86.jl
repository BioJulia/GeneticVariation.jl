# ng86.jl
# =======
#
# NG86 based distance measures for codons.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE.md

const CDN = Union{DNACodon, RNACodon}

"""
    classify_neighbor(codon::DNACodon)

Computes and classifies the neighbors of a given `codon` as either a
transition neighbor, or a transversion neighbor.
"""
function classify_neighbors(codon::DNACodon)
    tsn = Vector{DNACodon}()
    tvn = Vector{DNACodon}()
    codon_bases = collect(codon)
    @inbounds for n in 1:3
        i = codon_bases[n]
        for j in ACGT
            if i != j
                thiscdn = copy(codon_bases)
                thiscdn[n] = j
                ipur = ispurine(i)
                jpur = ispurine(j)
                topush = ifelse((ipur && jpur) || (!ipur && !jpur), tsn, tvn)
                push!(topush, DNACodon(thiscdn...))
            end
        end
    end
    return tsn, tvn
end

"""
    expected_ng86{C<:CDN}(codon::C, k::Float64 = 1.0, code::GeneticCode)

Enumerate the number of synonymous and non-synonymous sites present at a codon,
where each site may be both partially synonymous and non-synonymous.
"""
function expected_ng86{C<:CDN}(codon::C, k::Float64 = 1.0, code::BioSequences.GeneticCode)
    tsn, tvn = classify_neighbours(codon)
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
    for neighbor in tsv
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

function expected_ng86{C<:CDN}(codons::Vector{CDN}, k::Float64 = 1.0, code::BioSequences.GeneticCode)
    S = N = 0.0
    for codon in codons
        S_i, N_i = expected_ng86(NG86, codon, k, code)
        S += S_i
        N += N_i
    end
    return S, N
end

@inline function compare_codon{C<:CDN}(x::C, y::C, code::BioSequences.GeneticCode, weight::Float64 = 1.0)
    if code[x] == code[y]
        # Synonymous mutation.
        return weight, 0.0
    else
        # non-synonymous mutation.
        return 0.0, weight
    end
end

"""
    diff_positions{C<:CDN}(x::C, y::C)

Identify which sites in two codons are different.
"""
function diff_positions{C<:CDN}(x::C, y::C)
    diff_pos = Vector{Int}()
    for (i, k) in enumerate(zip(x, y))
        if k[1] != k[2]
            push!(diff_pos, i)
        end
    end
    return diff_pos, length(diff_pos)
end

function ng86_diff{C<:CDN}(x::C, y::C, code::GeneticCode)
    if x == y # Early escape, codons are the same, no syn or nonsyn mutations.
        return 0.0, 0.0
    else
        S = N = 0.0

        diff_pos, n_diffs = diff_positions(x, y) # Which positions are different.

        if n_diffs == 1

            # One site in the two codons is different. It is obvious and simple
            # then to count whether it is a synonymous or nonsynonymous mutation.
            S, N = compare_codon(x, y, code)
            return S, N

        elseif n_diffs == 2
            S = N = 0.0

            # For two changes, the number of synonymous and non-synonymous
            # differences per codon, sum to 2, there are two pathways,
            # each possible pathway having two steps.

            # For example, comparing CTA and GTT, the possible pathways are:
            # 1: CTA (L) -> GTA (V) -> GTT (V) : 1 nonsynonymous change and 1 synonymous change.
            # 2: CTA (L) -> CTT (L) -> GTT (V) : 1 nonsynonymous change and 1 synonymous change.

            @inbounds for pos in diff_pos
                bases = collect(x)
                bases[pos] = y[pos]
                temp_cdn = C(bases...)
                # Step 1 of pathway.
                S_i, N_i = compare_codon(x, temp_cdn, code, 0.5)
                S += S_i
                N += N_i
                # Step 2 of pathway.
                S_i, N_i = compare_codon(temp_cdn, y, code, 0.5)
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
                S_i, N_i = compare_codon(x, tmp_cdn_a, code, 0.5 / 3)
                S += S_i
                N += N_i
                S_i, N_i = compare_codon(tmp_cdn_a, tmp_cdn_b, code, 0.5 / 3)
                S += S_i
                N += N_i
                S_i, N_i = compare_codon(tmp_cdn_b, y, code, 0.5 / 3)
                S += S_i
                N += N_i
            end
        end

        return S, N
    end
end


function dNdS{C<:CDN}(x::Vector{C}, y::Vector{C}, k::Float64 = 1, code::GeneticCode)
    S_sites_x, N_sites_x = ng86_sites(x, k, code)
    S_sites_y, N_sites_y = ng86_sites(y, k, code)
    S_sites = (S_sites_x + S_sites_y) / 2
    N_sites = (N_sites_x + N_sites_y) / 2
    S = N = 0.0
    for (i, j) in zip(x, y)
        S_i, N_i = ng86_diff(i, j, code)

    end
end
