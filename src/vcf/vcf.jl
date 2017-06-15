# vcf.jl
# ===================
#
# A submodule for reading, writing, and working with VCF formatted files.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE

module VCF

import Automa
import Automa.RegExp: @re_str
import BioCore: BioCore, isfilled
import BioCore.Exceptions: missingerror
import BufferedStreams

include("record.jl")
include("metainfo.jl")
include("header.jl")
include("reader.jl")
include("writer.jl")

end
