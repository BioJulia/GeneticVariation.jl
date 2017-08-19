```@meta
CurrentModule = GeneticVariation
DocTestSetup = quote
    using GeneticVariation
end
```

# Computing allele frequencies

When first looking at the diversity present in a population, it is common to
want to know how many of each unique allele there is in a population i.e. the
allele frequencies of the population are.

Formally defined, allele frequency is a measure of the relative frequency of an
allele on a genetic locus in a population.

In population genetics, allele frequencies show the genetic diversity of a
species population or equivalently the richness of its gene pool.

Population genetics studies the different "forces" that might lead to changes
in the distribution and frequencies of alleles - in other words, to evolution.

Besides selection, these forces include genetic drift, mutation and migration.

Computing allele frequencies then, is an essential task for many wishing to
work with genetic variation, and so methods for computing such frequencies
are included in GeneticVariation.jl.

Allele frequencies can be computed for genes, micro-satellites, RFPL patterns,
and from SNPs.

```@docs
gene_frequencies
```

# Computing measures of genetic diversity

## Nucleotide diversity

```@docs
nuc_div
```
