```@meta
CurrentModule = GeneticVariation
DocTestSetup = quote
    using GeneticVariation
end
```

# Distances

GeneticVariation.jl allows one to compute a range of distance measures that
are relevant to the analysis of genetic variation.

## MASH distances

```@docs
MASH
```

You can generate a MinHash sketch using the `minhash()` function in our
`BiologicalSequences.jl` package.

```jlcon
using BiologicalSequences

seq1 = dna"AAATAAGGCACAACTATGCAT"
sketch1 = minhash(seq, 5, 10)
```

Then, if you have MinHash sketches with the same parameters for two sequences,
you can determine the MASH distance between them.

```jlcon
seq2 = dna"AATTAACGCACGGACTGCGGTAAT"
sketch2 = minhash(seq, 5, 10)

using GeneticVariation

distance(MASH, sketch1, sketch2)
# or use the convenience function
mash(sketch1, sketch2)
```

For more information on what size kmers and what size sketches are appropriate
for your use-case, see [Odnov et. al.](http://doi.org/10.1186/s13059-016-0997-x)
in _Genome Biology_.
