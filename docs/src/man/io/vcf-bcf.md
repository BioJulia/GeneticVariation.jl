```@meta
CurrentModule = GeneticVariation
DocTestSetup = quote
    using GeneticVariation
end
```

# VCF and BCF Formatted files

VCF is a text-based file format for representing genetic polymorphism.

VCF files can be read using `VCFReader`, respectively:

```julia
reader = VCF.Reader(open("example.vcf", "r"))
for record in reader
    # do something
end
close(reader)
```

A reader first reads the header section of a file and creates a `VCF.Header`
object. The `header` function is provided to access the header object of a
reader:

```jlcon
julia> header(reader)
GeneticVariation.VCF.Header:
  metainfo tags: fileformat fileDate source reference contig phasing INFO FILTER FORMAT
     sample IDs: NA00001 NA00002 NA00003

julia> find(header(reader), "FORMAT")
4-element Array{GeneticVariation.VCF.MetaInfo,1}:
 GeneticVariation.VCF.MetaInfo:
    tag: FORMAT
  value: ID="GT" Number="1" Type="String" Description="Genotype"          
 GeneticVariation.VCF.MetaInfo:
    tag: FORMAT
  value: ID="GQ" Number="1" Type="Integer" Description="Genotype Quality"
 GeneticVariation.VCF.MetaInfo:
    tag: FORMAT
  value: ID="DP" Number="1" Type="Integer" Description="Read Depth"       
 GeneticVariation.VCF.MetaInfo:
    tag: FORMAT
  value: ID="HQ" Number="2" Type="Integer" Description="Haplotype Quality"
```

`VCF.MetaInfo` variables in the header support the following accessors:

| Accessor      | Description                          |
| :-------      | :----------                          |
| `metainfotag` | tag string                           |
| `metainfoval` | value string                         |
| `keys`        | keys of fields between '<' and '>'   |
| `values`      | values of fields between '<' and '>' |
| `[<key>]`     | value of a field with `key`          |

```jlcon
julia> metainfo = VCF.MetaInfo("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
GeneticVariation.VCF.MetaInfo:
    tag: FORMAT
  value: ID="GT" Number="1" Type="String" Description="Genotype"

julia> metainfotag(metainfo)
"FORMAT"

julia> metainfoval(metainfo)
"<ID=GT,Number=1,Type=String,Description=\"Genotype\">"

julia> keys(metainfo)
4-element Array{String,1}:
 "ID"         
 "Number"     
 "Type"       
 "Description"

julia> metainfo["ID"]
"GT"

```

`VCF.Record` and `BCF.Record` variables support the following accessor functions
(see the docstring of each accessor for the details):

| Accessor   | Description                    |
| :-------   | :----------                    |
| `chrom`    | chromosome name                |
| `pos`      | reference position             |
| `id`       | unique identifiers             |
| `ref`      | reference bases                |
| `alt`      | alternate bases                |
| `qual`     | Phred-scaled quality score     |
| `filter`   | filter status                  |
| `info`     | additional information         |
| `infokeys` | keys of additional information |
| `format`   | genotype format                |
| `genotype` | genotype information           |

```jlcon
julia> record = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
GeneticVariation.VCF.Record:
   chromosome: 20
     position: 14370
   identifier: rs6054257
    reference: G
    alternate: A
      quality: 29.0
       filter: PASS
  information: NS=3 DP=14 AF=0.5 DB H2
       format: GT GQ DP HQ
     genotype: [1] 0|0 48 1 51,51 [2] 1|0 48 8 51,51

julia> VCF.chrom(record)
"20"

julia> VCF.pos(record)
14370

julia> VCF.id(record)
1-element Array{String,1}:
 "rs6054257"

julia> VCF.ref(record)
"G"

julia> VCF.alt(record)
1-element Array{String,1}:
 "A"

julia> VCF.qual(record)
29.0

julia> VCF.filter(record)
1-element Array{String,1}:
 "PASS"

julia> VCF.info(record)
5-element Array{Pair{String,String},1}:
 "NS"=>"3"  
 "DP"=>"14"
 "AF"=>"0.5"
 "DB"=>""   
 "H2"=>""   

julia> VCF.format(record)
4-element Array{String,1}:
 "GT"
 "GQ"
 "DP"
 "HQ"

julia> VCF.genotype(record)
2-element Array{Array{String,1},1}:
 String["0|0","48","1","51,51"]
 String["1|0","48","8","51,51"]

julia> VCF.genotype(record, 1:2, "GT")
2-element Array{String,1}:
 "0|0"
 "1|0"

julia> VCF.genotype(record, 1:1, "GT")
1-element Array{String,1}:
 "0|0"

julia> VCF.genotype(record, 1:2, "GT")
2-element Array{String,1}:
 "0|0"
 "1|0"

```
