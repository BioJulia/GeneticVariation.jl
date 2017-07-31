# writer.jl
# =========
#
# A writer for BCF formatted files.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE

struct Writer{T<:IO} <: BioCore.IO.AbstractWriter
    stream::BGZFStreams.BGZFStream{T}
end

"""
    BCF.Writer(output::IO, header::VCF.Header)
Create a data writer of the BCF file format.
# Arguments
* `output`: data sink
* `header`: VCF header object
"""
function Writer(output::IO, header::VCF.Header)
    stream = BGZFStreams.BGZFStream(output, "w")
    write(stream, b"BCF\x02\x02")
    buf = IOBuffer()
    len = write(buf, header)
    if len > typemax(Int32)
        error("too long header")
    end
    write(stream, htol(Int32(len)))
    data = take!(buf)
    @assert length(data) == len
    write(stream, data)
    return Writer(stream)
end

function BioCore.IO.stream(writer::Writer)
    return writer.stream
end

function Base.write(writer::Writer, record::Record)
    n = 0
    n += write(writer.stream, htol(record.sharedlen))
    n += write(writer.stream, htol(record.indivlen))
    n += write(writer.stream, record.data)
    return n
end
