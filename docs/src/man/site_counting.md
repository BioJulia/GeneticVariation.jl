```@meta
CurrentModule = GeneticVariation
DocTestSetup = quote
    using GeneticVariation
end
```

# Identifying and counting sequence sites

GeneticVariation.jl extends the [site-counting](https://biojulia.github.io/BioSequences.jl/stable/sequences/bioseq/#site-counting)
methods in BioSequences.jl, using the same fast bit-parallel techniques to rapidly
compute the numbers of different types of mutations between two large biological sequences.
Such computation is required for many population genetic analyses of variation, such
as computation of evolutionary distances.

## Types of site added

```@docs
Conserved
Mutated
Segregating
```

See the [site-counting](https://biojulia.github.io/BioSequences.jl/stable/sequences/bioseq/#site-counting)
section of the BioSequences.jl documentation to see how to use the `count` and
`count_pairwise` methods to count different types of site.
