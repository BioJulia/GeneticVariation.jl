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
    distance,
    Jaccard,
    MASH,
    mash,
    jaccard,

    # VCF and BCF
    VCF,
    BCF,
    header,
    metainfotag,
    metainfoval,
    isfilled,
    MissingFieldException
    #mashdistance,
    #distance,
    #Proportion

importall BioCore
import BioSymbols: ispurine, ispyrimidine
import BioSequences:
    BioSequences,
    Alphabet,
    DNAAlphabet,
    RNAAlphabet,
    BioSequence,
    MinHashSketch,
    Certain,
    Mismatch,
    Match,
    Position,
    bp_chunk_count
import Compat: @compat
import IntervalTrees: Interval, IntervalValue
import Twiddle:
    enumerate_nibbles,
    nibble_mask,
    count_zero_nibbles,
    count_nonzero_nibbles,
    count_one_nibbles,
    count_zero_bitpairs,
    count_nonzero_bitpairs

#include("distances/dist.jl")
include("vcf/vcf.jl")
include("bcf/bcf.jl")
include("site_counting.jl")
include("distances/minhash.jl")

end # Module GeneticVariation
