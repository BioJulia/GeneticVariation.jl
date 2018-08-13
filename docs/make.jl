using Documenter, GeneticVariation

makedocs(
    modules = [GeneticVariation],
    format = :html,
    sitename = "GeneticVariation.jl",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "I/O" => [
                "VCF and BCF formatted files" => "man/io/vcf-bcf.md"
            ],
            "Mutation Counting" => "man/site_counting.md",
            "Genetic Diversity" => "man/diversity.md"
        ]
    ],
    authors = "Kenta Sato, Ben J. Ward, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/GeneticVariation.jl.git",
    julia = "1.0",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing
)
