@testset "VCF" begin
    metainfo = VCF.MetaInfo()
    @test !isfilled(metainfo)
    @test ismatch(r"^GeneticVariation.VCF.MetaInfo: <not filled>", repr(metainfo))
    @test_throws ArgumentError metainfotag(metainfo)

    metainfo = VCF.MetaInfo(b"##source=foobar1234")
    @test isfilled(metainfo)
    @test metainfotag(metainfo) == "source"
    @test metainfoval(metainfo) == "foobar1234"

    metainfo = VCF.MetaInfo("##source=foobar1234")
    @test isfilled(metainfo)
    @test metainfotag(metainfo) == "source"
    @test metainfoval(metainfo) == "foobar1234"

    metainfo = VCF.MetaInfo(metainfo)
    @test isa(metainfo, VCF.MetaInfo)
    metainfo = VCF.MetaInfo(metainfo, tag="date")
    @test metainfotag(metainfo) == "date"
    metainfo = VCF.MetaInfo(metainfo, value="2017-01-30")
    @test metainfoval(metainfo) == "2017-01-30"
    metainfo = VCF.MetaInfo(metainfo, tag="INFO", value=["ID"=>"DP", "Number"=>"1", "Type"=>"Integer", "Description"=>"Total Depth"])
    @test metainfo["ID"] == "DP"
    @test metainfo["Number"] == "1"
    @test metainfo["Type"] == "Integer"
    @test metainfo["Description"] == "Total Depth"
    @test metainfotag(metainfo) == "INFO"
    @test metainfoval(metainfo) == """<ID=DP,Number=1,Type=Integer,Description="Total Depth">"""

    record = VCF.Record()
    @test !isfilled(record)
    @test ismatch(r"^GeneticVariation.VCF.Record: <not filled>", repr(record))
    @test_throws ArgumentError VCF.chrom(record)

    record = VCF.Record("20\t302\t.\tT\tTA\t999\t.\t.\tGT")
    @test isfilled(record)
    @test VCF.haschrom(record)
    @test VCF.chrom(record) == "20"
    @test VCF.haspos(record)
    @test VCF.pos(record) == 302
    @test !VCF.hasid(record)
    @test_throws MissingFieldException VCF.id(record)
    @test VCF.hasref(record)
    @test VCF.ref(record) == "T"
    @test VCF.hasalt(record)
    @test VCF.alt(record) == ["TA"]
    @test VCF.hasqual(record)
    @test VCF.qual(record) == 999
    @test !VCF.hasfilter(record)
    @test_throws MissingFieldException VCF.filter(record)
    @test VCF.infokeys(record) == String[]
    @test !VCF.hasinfo(record)
    @test_throws MissingFieldException VCF.info(record)
    @test VCF.hasformat(record)
    @test VCF.format(record) == ["GT"]

    # empty data is not a valid VCF record
    @test_throws ArgumentError VCF.Record("")
    @test_throws ArgumentError VCF.Record(b"")

    record = VCF.Record(b".\t.\t.\t.\t.\t.\t.\t.\t")
    @test isfilled(record)
    @test !VCF.haschrom(record)
    @test !VCF.haspos(record)
    @test !VCF.hasid(record)
    @test !VCF.hasref(record)
    @test !VCF.hasalt(record)
    @test !VCF.hasqual(record)
    @test !VCF.hasfilter(record)
    @test !VCF.hasinfo(record)

    record = VCF.Record(record)
    @test isa(record, VCF.Record)
    record = VCF.Record(record, chrom="chr1")
    @test VCF.chrom(record) == "chr1"
    record = VCF.Record(record, pos=1234)
    @test VCF.pos(record) == 1234
    record = VCF.Record(record, id="rs1111")
    @test VCF.id(record) == ["rs1111"]
    record = VCF.Record(record, ref="A")
    @test VCF.ref(record) == "A"
    record = VCF.Record(record, alt=["AT"])
    @test VCF.alt(record) == ["AT"]
    record = VCF.Record(record, qual=11.2)
    @test VCF.qual(record) == 11.2
    record = VCF.Record(record, filter="PASS")
    @test VCF.filter(record) == ["PASS"]
    record = VCF.Record(record, info=Dict("DP" => 20, "AA" => "AT", "DB"=>nothing))
    @test VCF.info(record, "DP") == "20"
    @test VCF.info(record, "AA") == "AT"
    @test VCF.info(record, "DB") == ""
    @test VCF.infokeys(record) == ["DP", "AA", "DB"]
    record = VCF.Record(record, genotype=[Dict("GT" => "0/0", "DP" => [10,20])])
    @test VCF.format(record) == ["DP", "GT"]
    @test VCF.genotype(record) == [["10,20", "0/0"]]

    let header = VCF.Header()
        @test isempty(header)
        push!(header, "##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta")
        @test !isempty(header)
        @test length(header) == 1
        unshift!(header, "##fileformat=VCFv4.3")
        @test length(header) == 2
        @test collect(header) == [
            VCF.MetaInfo("##fileformat=VCFv4.3"),
            VCF.MetaInfo("##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta")]
        @test startswith(repr(header), "GeneticVariation.VCF.Header:")
    end

    let header = VCF.Header(["##fileformat=VCFv4.3"], ["Sample1"])
        @test !isempty(header)
        @test length(header) == 1
        @test header.sampleID == ["Sample1"]
        @test first(header) == VCF.MetaInfo("##fileformat=VCFv4.3")
    end

    # minimum header
    data = b"""
    ##fileformat=VCFv4.3
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    """
    reader = VCF.Reader(BufferedInputStream(data))
    @test isa(header(reader), VCF.Header)
    let header = header(reader)
        @test length(header.metainfo) == 1
        @test metainfotag(header.metainfo[1]) == "fileformat"
        @test metainfoval(header.metainfo[1]) == "VCFv4.3"
        @test isempty(header.sampleID)
    end

    # realistic header
    data = b"""
    ##fileformat=VCFv4.2
    ##fileDate=20090805
    ##source=myImputationProgramV3.1
    ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
    ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
    ##phasing=partial
    ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
    """
    reader = VCF.Reader(BufferedInputStream(data))
    @test isa(header(reader), VCF.Header)

    let header = header(reader)
        @test length(header.metainfo) == 10

        let metainfo = header.metainfo[1]
            @test metainfotag(metainfo) == "fileformat"
            @test metainfoval(metainfo) == "VCFv4.2"
            @test_throws ArgumentError keys(metainfo)
            @test_throws ArgumentError values(metainfo)
        end
        @test length(find(header, "fileformat")) == 1
        @test first(find(header, "fileformat")) == header.metainfo[1]

        let metainfo = header.metainfo[2]
            @test metainfotag(metainfo) == "fileDate"
            @test metainfoval(metainfo) == "20090805"
            @test_throws ArgumentError keys(metainfo)
            @test_throws ArgumentError values(metainfo)
        end

        let metainfo = header.metainfo[5]
            @test metainfotag(metainfo) == "contig"
            @test metainfoval(metainfo) == """<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>"""
            @test keys(metainfo) == ["ID", "length", "assembly", "md5", "species", "taxonomy"]
            @test values(metainfo) == ["20", "62435964", "B36", "f126cdf8a6e0c7f379d618ff66beb2da", "Homo sapiens", "x"]
            @test metainfo["ID"] == "20"
            @test metainfo["md5"] == "f126cdf8a6e0c7f379d618ff66beb2da"
            @test metainfo["taxonomy"] == "x"
        end

        let metainfo = header.metainfo[7]
            @test metainfotag(metainfo) == "INFO"
            @test metainfoval(metainfo) == """<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">"""
            @test keys(metainfo) == ["ID", "Number", "Type", "Description"]
            @test values(metainfo) == ["NS", "1", "Integer", "Number of Samples With Data"]
            @test metainfo["ID"] == "NS"
            @test metainfo["Type"] == "Integer"
        end
        @test length(find(header, "INFO")) == 4

        @test header.sampleID == ["NA00001", "NA00002", "NA00003"]
    end

    data = b"""
    ##fileformat=VCFv4.3
    ##contig=<ID=chr1>
    ##contig=<ID=chr2>
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    ##FORMAT=<ID=GT,Number=1,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002
    chr1\t1234\trs001234\tA\tC\t30\tPASS\tDP=10;AF=0.3\tGT\t0|0\t0/1
    chr2\t4\t.\tA\tAA,AAT\t.\t.\tDP=5\tGT:DP\t0|1:42\t0/1
    """
    reader = VCF.Reader(BufferedInputStream(data))
    record = VCF.Record()

    @test read!(reader, record) === record
    @test VCF.chrom(record) == "chr1"
    @test VCF.pos(record) === 1234
    @test VCF.id(record) == ["rs001234"]
    @test VCF.ref(record) == "A"
    @test VCF.alt(record) == ["C"]
    @test VCF.qual(record) === 30.0
    @test VCF.filter(record) == ["PASS"]
    @test VCF.info(record) == ["DP" => "10", "AF" => "0.3"]
    @test VCF.info(record, "DP") == "10"
    @test VCF.info(record, "AF") == "0.3"
    @test VCF.format(record) == ["GT"]
    @test VCF.genotype(record) == [["0|0"], ["0/1"]]
    @test VCF.genotype(record, 1) == VCF.genotype(record)[1]
    @test VCF.genotype(record, 2) == VCF.genotype(record)[2]
    @test VCF.genotype(record, 1, "GT") == "0|0"
    @test VCF.genotype(record, 2, "GT") == "0/1"
    @test VCF.genotype(record, 1:2, "GT") == ["0|0", "0/1"]
    @test VCF.genotype(record, :, "GT") == VCF.genotype(record, 1:2, "GT")
    @test ismatch(r"^GeneticVariation.VCF.Record:\n.*", repr(record))

    @test read!(reader, record) === record
    @test VCF.chrom(record) == "chr2"
    @test VCF.pos(record) == 4
    @test !VCF.hasid(record)
    @test VCF.ref(record) == "A"
    @test VCF.alt(record) == ["AA", "AAT"]
    @test !VCF.hasqual(record)
    @test !VCF.hasfilter(record)
    @test VCF.info(record) == ["DP" => "5"]
    @test VCF.info(record, "DP") == "5"
    @test_throws KeyError VCF.info(record, "AF")
    @test VCF.format(record) == ["GT", "DP"]
    @test VCF.genotype(record) == [["0|1", "42"], ["0/1", "."]]
    @test VCF.genotype(record, 1) == VCF.genotype(record)[1]
    @test VCF.genotype(record, 2) == VCF.genotype(record)[2]
    @test VCF.genotype(record, 1, "GT") == "0|1"
    @test VCF.genotype(record, 1, "DP") == "42"
    @test VCF.genotype(record, 2, "GT") == "0/1"
    @test VCF.genotype(record, 2, "DP") == "."
    @test VCF.genotype(record, 1:2, "GT") == ["0|1", "0/1"]
    @test VCF.genotype(record, 1:2, "DP") == ["42", "."]
    @test VCF.genotype(record, :, "DP") == VCF.genotype(record, 1:2, "DP")
    @test_throws KeyError VCF.genotype(record, :, "BAD")

    @test_throws EOFError read!(reader, record)

    # round-trip test
    vcfdir = joinpath(fmtdir, "VCF")
    for specimen in YAML.load_file(joinpath(vcfdir, "index.yml"))
        filepath = joinpath(vcfdir, specimen["filename"])
        records = VCF.Record[]
        reader = open(VCF.Reader, filepath)
        output = IOBuffer()
        writer = VCF.Writer(output, header(reader))
        for record in reader
            write(writer, record)
            push!(records, record)
        end
        close(reader)
        flush(writer)

        records2 = VCF.Record[]
        for record in VCF.Reader(IOBuffer(takebuf_array(output)))
            push!(records2, record)
        end
        @test records == records2
    end
end

function parsehex(str)
    return map(x -> parse(UInt8, x, 16), split(str, ' '))
end
