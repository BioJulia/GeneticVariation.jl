# GeneticVariation.jl

[![Latest Release](https://img.shields.io/github/release/BioJulia/GeneticVariation.jl.svg)](https://github.com/BioJulia/GeneticVariation.jl/releases/latest)
[![GeneticVariation](http://pkg.julialang.org/badges/GeneticVariation_0.6.svg)](http://pkg.julialang.org/?pkg=GeneticVariation)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE)
![BioJulia maintainer: bicycle1885](https://img.shields.io/badge/BioJulia%20Maintainer-bicycle1885-orange.svg)
![BioJulia maintainer: Ward9250](https://img.shields.io/badge/BioJulia%20Maintainer-Ward9250-orange.svg)

**Development builds:**
[![Build Status](https://travis-ci.org/BioJulia/GeneticVariation.jl.svg?branch=master)](https://travis-ci.org/BioJulia/GeneticVariation.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/29um8ekg6en3s23a?svg=true)](https://ci.appveyor.com/project/Ward9250/geneticvariation-jl)
[![codecov](https://codecov.io/gh/BioJulia/GeneticVariation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/BioJulia/GeneticVariation.jl)

## Description

GeneticVariation.jl provides types and methods for working with genetic variation.
It provides a VCF and BCF parser, as well as methods for working with variation
in sequences such as evolutionary distance computation, and counting different
mutation types.

## Installation

Install GeneticVariation from the Julia REPL:

```julia
julia> Pkg.add("GeneticVariation")
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.
