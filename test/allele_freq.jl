@testset "Allele frequenies" begin
    @testset "Gene frequencies" begin

        sequences = ["ATCGATCG", "AGGGGG", "CCCCCCCCCCCCCCC", "TTTTTCCCC",
                     "ATCGATCG", "AGGGGG", "ATCGATCG", "ATCGATCG"]

        answer = Dict{String, Float64}("TTTTTCCCC" => 0.125,
                                       "CCCCCCCCCCCCCCC" => 0.125,
                                       "ATCGATCG" => 0.5,
                                       "AGGGGG" => 0.25)

        function test_gene_frequencies(genes, answer, seqtype)
            test_genes = convert(Vector{seqtype}, sequences)
            test_answer = convert(Dict{seqtype, Float64}, answer)
            @test gene_frequencies(test_genes) == test_answer
            @test GeneticVariation._gene_frequencies(test_genes, seqtype, Base.SizeUnknown()) == test_answer
            @test GeneticVariation._gene_frequencies(test_genes, seqtype, Base.HasShape()) == test_answer
        end

        for n in (2, 4)
            test_gene_frequencies(sequences, answer, BioSequence{DNAAlphabet{n}})
        end
    end
end
