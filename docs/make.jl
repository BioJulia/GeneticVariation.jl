using Pkg, Documenter
using GeneticVariation

makedocs(
    checkdocs = :all,
    linkcheck = true,
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        edit_link = "master",
        sidebar_sitename = false
    ),
    modules = [GeneticVariation],
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
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ) * ", The BioJulia Organisation, and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/GeneticVariation.jl.git",
    devbranch = "master",
    push_preview = true
)
