# bcf.jl
# ======
#
# A submodule for reading, writing, and working with BCF formatted files.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE

module BCF

import BioCore: BioCore, isfilled
import BGZFStreams
import BufferedStreams

include("record.jl")
include("reader.jl")
include("writer.jl")

end
