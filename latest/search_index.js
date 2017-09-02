var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#GeneticVariation.jl-1",
    "page": "Home",
    "title": "GeneticVariation.jl",
    "category": "section",
    "text": "(Image: Latest Release) (Image: GeneticVariation) (Image: GeneticVariation) (Image: License) (Image: BioJulia maintainer: bicycle1885) (Image: BioJulia maintainer: Ward9250)Development builds: (Image: Build Status) (Image: Build status) (Image: codecov)"
},

{
    "location": "index.html#Description-1",
    "page": "Home",
    "title": "Description",
    "category": "section",
    "text": "GeneticVariation.jl provides types and methods for working with genetic variation. It provides a VCF and BCF parser, as well as methods for working with variation in sequences such as evolutionary distance computation, and counting different mutation types."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Install GeneticVariation from the Julia REPL:julia> Pkg.add(\"GeneticVariation\")If you are interested in the cutting edge of the development, please check out the master branch to try new features before release."
},

{
    "location": "man/io/vcf-bcf.html#",
    "page": "VCF and BCF formatted files",
    "title": "VCF and BCF formatted files",
    "category": "page",
    "text": "CurrentModule = GeneticVariation\nDocTestSetup = quote\n    using GeneticVariation\nend"
},

{
    "location": "man/io/vcf-bcf.html#VCF-and-BCF-Formatted-files-1",
    "page": "VCF and BCF formatted files",
    "title": "VCF and BCF Formatted files",
    "category": "section",
    "text": "VCF is a text-based file format for representing genetic polymorphism.VCF files can be read using VCFReader, respectively:reader = VCF.Reader(open(\"example.vcf\", \"r\"))\nfor record in reader\n    # do something\nend\nclose(reader)A reader first reads the header section of a file and creates a VCF.Header object. The header function is provided to access the header object of a reader:julia> header(reader)\nGeneticVariation.VCF.Header:\n  metainfo tags: fileformat fileDate source reference contig phasing INFO FILTER FORMAT\n     sample IDs: NA00001 NA00002 NA00003\n\njulia> find(header(reader), \"FORMAT\")\n4-element Array{GeneticVariation.VCF.MetaInfo,1}:\n GeneticVariation.VCF.MetaInfo:\n    tag: FORMAT\n  value: ID=\"GT\" Number=\"1\" Type=\"String\" Description=\"Genotype\"          \n GeneticVariation.VCF.MetaInfo:\n    tag: FORMAT\n  value: ID=\"GQ\" Number=\"1\" Type=\"Integer\" Description=\"Genotype Quality\"\n GeneticVariation.VCF.MetaInfo:\n    tag: FORMAT\n  value: ID=\"DP\" Number=\"1\" Type=\"Integer\" Description=\"Read Depth\"       \n GeneticVariation.VCF.MetaInfo:\n    tag: FORMAT\n  value: ID=\"HQ\" Number=\"2\" Type=\"Integer\" Description=\"Haplotype Quality\"VCF.MetaInfo variables in the header support the following accessors:Accessor Description\nmetainfotag tag string\nmetainfoval value string\nkeys keys of fields between '<' and '>'\nvalues values of fields between '<' and '>'\n[<key>] value of a field with keyjulia> metainfo = VCF.MetaInfo(\"##FORMAT=<ID=GT,Number=1,Type=String,Description=\\\"Genotype\\\">\")\nGeneticVariation.VCF.MetaInfo:\n    tag: FORMAT\n  value: ID=\"GT\" Number=\"1\" Type=\"String\" Description=\"Genotype\"\n\njulia> metainfotag(metainfo)\n\"FORMAT\"\n\njulia> metainfoval(metainfo)\n\"<ID=GT,Number=1,Type=String,Description=\\\"Genotype\\\">\"\n\njulia> keys(metainfo)\n4-element Array{String,1}:\n \"ID\"         \n \"Number\"     \n \"Type\"       \n \"Description\"\n\njulia> metainfo[\"ID\"]\n\"GT\"\nVCF.Record and BCF.Record variables support the following accessor functions (see the docstring of each accessor for the details):Accessor Description\nchrom chromosome name\npos reference position\nid unique identifiers\nref reference bases\nalt alternate bases\nqual Phred-scaled quality score\nfilter filter status\ninfo additional information\ninfokeys keys of additional information\nformat genotype format\ngenotype genotype informationjulia> record = VCF.Record(\"20\\t14370\\trs6054257\\tG\\tA\\t29\\tPASS\\tNS=3;DP=14;AF=0.5;DB;H2\\tGT:GQ:DP:HQ\\t0|0:48:1:51,51\\t1|0:48:8:51,51\")\nGeneticVariation.VCF.Record:\n   chromosome: 20\n     position: 14370\n   identifier: rs6054257\n    reference: G\n    alternate: A\n      quality: 29.0\n       filter: PASS\n  information: NS=3 DP=14 AF=0.5 DB H2\n       format: GT GQ DP HQ\n     genotype: [1] 0|0 48 1 51,51 [2] 1|0 48 8 51,51\n\njulia> VCF.chrom(record)\n\"20\"\n\njulia> VCF.pos(record)\n14370\n\njulia> VCF.id(record)\n1-element Array{String,1}:\n \"rs6054257\"\n\njulia> VCF.ref(record)\n\"G\"\n\njulia> VCF.alt(record)\n1-element Array{String,1}:\n \"A\"\n\njulia> VCF.qual(record)\n29.0\n\njulia> VCF.filter(record)\n1-element Array{String,1}:\n \"PASS\"\n\njulia> VCF.info(record)\n5-element Array{Pair{String,String},1}:\n \"NS\"=>\"3\"  \n \"DP\"=>\"14\"\n \"AF\"=>\"0.5\"\n \"DB\"=>\"\"   \n \"H2\"=>\"\"   \n\njulia> VCF.format(record)\n4-element Array{String,1}:\n \"GT\"\n \"GQ\"\n \"DP\"\n \"HQ\"\n\njulia> VCF.genotype(record)\n2-element Array{Array{String,1},1}:\n String[\"0|0\",\"48\",\"1\",\"51,51\"]\n String[\"1|0\",\"48\",\"8\",\"51,51\"]\n\njulia> VCF.genotype(record, 1:2, \"GT\")\n2-element Array{String,1}:\n \"0|0\"\n \"1|0\"\n\njulia> VCF.genotype(record, 1:1, \"GT\")\n1-element Array{String,1}:\n \"0|0\"\n\njulia> VCF.genotype(record, 1:2, \"GT\")\n2-element Array{String,1}:\n \"0|0\"\n \"1|0\"\n"
},

{
    "location": "man/site_counting.html#",
    "page": "Mutation Counting",
    "title": "Mutation Counting",
    "category": "page",
    "text": "CurrentModule = GeneticVariation\nDocTestSetup = quote\n    using GeneticVariation\nend"
},

{
    "location": "man/site_counting.html#Identifying-and-counting-sequence-sites-1",
    "page": "Mutation Counting",
    "title": "Identifying and counting sequence sites",
    "category": "section",
    "text": "GeneticVariation.jl extends the site-counting methods in BioSequences.jl, using the same fast bit-parallel techniques to rapidly compute the numbers of different types of mutations between two large biological sequences. Such computation is required for many population genetic analyses of variation, such as computation of evolutionary distances."
},

{
    "location": "man/site_counting.html#GeneticVariation.Conserved",
    "page": "Mutation Counting",
    "title": "GeneticVariation.Conserved",
    "category": "Type",
    "text": "A Conserved site describes a site where two aligned nucleotides are definately conserved. By definately conserved this means that the symbols of the site are non-ambiguity symbols, and they are the same symbol.\n\n\n\n"
},

{
    "location": "man/site_counting.html#GeneticVariation.Mutated",
    "page": "Mutation Counting",
    "title": "GeneticVariation.Mutated",
    "category": "Type",
    "text": "A Mutated site describes a site where two aligned nucleotides are definately mutated. By definately mutated this means that the symbols of the site are non-ambiguity symbols, and they are not the same symbol.\n\n\n\n"
},

{
    "location": "man/site_counting.html#Types-of-site-added-1",
    "page": "Mutation Counting",
    "title": "Types of site added",
    "category": "section",
    "text": "Conserved\nMutatedSee the [site-counting](site-counting section of the BioSequences.jl documentation to see how to use the count and count_pairwise methods to count different types of site."
},

{
    "location": "man/diversity.html#",
    "page": "Genetic Diversity",
    "title": "Genetic Diversity",
    "category": "page",
    "text": "CurrentModule = GeneticVariation\nDocTestSetup = quote\n    using GeneticVariation\nend"
},

{
    "location": "man/diversity.html#GeneticVariation.gene_frequencies",
    "page": "Genetic Diversity",
    "title": "GeneticVariation.gene_frequencies",
    "category": "Function",
    "text": "gene_frequencies(seqcounts::Composition{T}) where T <: Sequence\n\nCompute gene_frequencies from a BioSequences.Composition variable that contains unique sequence counts.\n\n\n\ngene_frequencies(iterable)\n\nCompute the gene frequencies for any iterable with an eltype which is a concrete subtype of the abstract Sequence type.\n\n\n\n"
},

{
    "location": "man/diversity.html#Computing-allele-frequencies-1",
    "page": "Genetic Diversity",
    "title": "Computing allele frequencies",
    "category": "section",
    "text": "When first looking at the diversity present in a population, it is common to want to know how many of each unique allele there is in a population i.e. the allele frequencies of the population are.Formally defined, allele frequency is a measure of the relative frequency of an allele on a genetic locus in a population.In population genetics, allele frequencies show the genetic diversity of a species population or equivalently the richness of its gene pool.Population genetics studies the different \"forces\" that might lead to changes in the distribution and frequencies of alleles - in other words, to evolution.Besides selection, these forces include genetic drift, mutation and migration.Computing allele frequencies then, is an essential task for many wishing to work with genetic variation, and so methods for computing such frequencies are included in GeneticVariation.jl.Allele frequencies can be computed for genes, micro-satellites, RFPL patterns, and from SNPs.gene_frequencies"
},

{
    "location": "man/diversity.html#Computing-measures-of-genetic-diversity-1",
    "page": "Genetic Diversity",
    "title": "Computing measures of genetic diversity",
    "category": "section",
    "text": ""
},

{
    "location": "man/diversity.html#GeneticVariation.NL79",
    "page": "Genetic Diversity",
    "title": "GeneticVariation.NL79",
    "category": "Function",
    "text": "NL79(m::M, f::V) where {M<:AbstractMatrix{Float64},V<:AbstractVector{Float64}}\n\nCompute nucleotide diversity using a matrix of the number of mutations between sequence pairs, and a vector of the frequencies of each sequence in the population.\n\n\n\nNL79(sequences)\n\nCompute nucleotide diversity, as described by Nei and Li in 1979.\n\nThis measure is defined as the average number of nucleotide differences per site between two DNA sequences in all possible pairs in the sample population, and is often denoted by the greek letter pi.\n\nSequences should be any iterable that yields biosequence types.\n\nExamples\n\njulia> testSeqs = [dna\"AAAACTTTTACCCCCGGGGG\",\n                   dna\"AAAACTTTTACCCCCGGGGG\",\n                   dna\"AAAACTTTTACCCCCGGGGG\",\n                   dna\"AAAACTTTTACCCCCGGGGG\",\n                   dna\"AAAAATTTTACCCCCGTGGG\",\n                   dna\"AAAAATTTTACCCCCGTGGG\",\n                   dna\"AAAACTTTTTCCCCCGTAGG\",\n                   dna\"AAAACTTTTTCCCCCGTAGG\",\n                   dna\"AAAAATTTTTCCCCCGGAGG\",\n                   dna\"AAAAATTTTTCCCCCGGAGG\"]\n10-element Array{BioSequences.BioSequence{BioSequences.DNAAlphabet{4}},1}:\n AAAACTTTTACCCCCGGGGG\n AAAACTTTTACCCCCGGGGG\n AAAACTTTTACCCCCGGGGG\n AAAACTTTTACCCCCGGGGG\n AAAAATTTTACCCCCGTGGG\n AAAAATTTTACCCCCGTGGG\n AAAACTTTTTCCCCCGTAGG\n AAAACTTTTTCCCCCGTAGG\n AAAAATTTTTCCCCCGGAGG\n AAAAATTTTTCCCCCGGAGG\n\n julia> NL79(testSeqs)\n 0.096\n\n\n\n\n"
},

{
    "location": "man/diversity.html#Nucleotide-diversity-1",
    "page": "Genetic Diversity",
    "title": "Nucleotide diversity",
    "category": "section",
    "text": "Nucleotide diversity is a concept in molecular genetics which is used to measure the degree of polymorphism within a population.There are different methods which can be used to compute measures of nucleotide diversity, we list them below, and show how to compute them using GeneticVariation. NL79"
},

{
    "location": "contributing.html#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "contributing.html#Contributing-1",
    "page": "Contributing",
    "title": "Contributing",
    "category": "section",
    "text": "We appreciate contributions from users including reporting bugs, fixing issues, improving performance and adding new features.If you have a question about contributing or using this package, you are encouraged to use the Bio category of the Julia discourse site.Detailed guidance for contributing to all BioJulia packages is provided at the BioJulia Contribution Documentation.Here we list specific details about contributing and maintainership pertaining specifically to the GeneticVariation.jl package."
},

{
    "location": "contributing.html#Named-maintainers-1",
    "page": "Contributing",
    "title": "Named maintainers",
    "category": "section",
    "text": "The named maintainers of this package are Kenta Sato and Ben Ward. It is their responsibility to make final choices about pull requests and issues, although because of our community structure, you will find other maintainers assisting them."
},

{
    "location": "contributing.html#Branching-model-1",
    "page": "Contributing",
    "title": "Branching model",
    "category": "section",
    "text": "The branching model used to develop and make releases of this package is the OneFlow model summarized in the BioJulia Contribution Documentation"
},

]}
