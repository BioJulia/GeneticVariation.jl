# GeneticVariation.jl
# ===================
#
# A julia package for the representation and analysis of genetic variation.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE

__precompile__()

module GeneticVariation

export
    # Site types
    Conserved,
    Mutated,
    #Transition,
    #Transversion,

    # Distances
    Proportion,
    Jaccard,
    MASH,
    distance,
    pdistance,
    mash,
    jaccard,
    NG86,

    # Allele frequencies
    gene_frequencies,

    # Nucleotide diversity
    NL79,

    # VCF and BCF
    VCF,
    BCF,
    header,
    metainfotag,
    metainfoval,
    isfilled,
    MissingFieldException

importall BioCore
import BioSequences:
    BioSequences,
    Alphabet,
    AA_Term,
    BioSequence,
    bp_chunk_count,
    Certain,
    Composition,
    DNAAlphabet,
    GeneticCode,
    ispurine,
    Kmer,
    Match,
    Mismatch,
    MinHashSketch,
    NucAlphs,
    Position,
    RNAAlphabet,
    Sequence

import Compat: @compat
import Combinatorics.permutations
import IntervalTrees: Interval, IntervalValue
import Twiddle:
    enumerate_nibbles,
    nibble_mask,
    count_zero_nibbles,
    count_nonzero_nibbles,
    count_one_nibbles,
    count_zero_bitpairs,
    count_nonzero_bitpairs

include("vcf/vcf.jl")
include("bcf/bcf.jl")
include("site_counting.jl")
include("distances/minhash.jl")
include("distances/proportion.jl")
include("dnds.jl")
include("allele_freq.jl")
include("nuc_div.jl")

end # Module GeneticVariation
