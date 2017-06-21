using Documenter, GeneticVariation

makedocs()
deploydocs(
    deps = Deps.pip("mkdocs", "pygments", "mkdocs-material"),
    repo = "github.com/BioJulia/GeneticVariation.jl.git",
    julia = "0.5",
    osname = "linux",
)
