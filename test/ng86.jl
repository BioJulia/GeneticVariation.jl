@testset "NG86" begin
    codonsA = [kmer"ATG",
               kmer"AAA",
               kmer"CCC",
               kmer"GGG",
               kmer"TTT"]

    codonsB = [kmer"ATG",
               kmer"AAA",
               kmer"CGC",
               kmer"GGC",
               kmer"TAC"]

    @testset "expected" begin
        n_ans = [3.0, 2.666, 2.0, 2.0, 2.666]
        s_ans = [0.0, 0.333, 1.0, 1.0, 0.333]
        for i in 1:endof(codonsA)
            cdn = codonsA[i]
            @test GeneticVariation.S_N_NG86(cdn, 1.0, ncbi_trans_table[1])[1] ≈ s_ans[i] atol=0.001
            @test GeneticVariation.S_N_NG86(cdn, 1.0, ncbi_trans_table[1])[2] ≈ n_ans[i] atol=0.001
        end
    end

    @testset "observed" begin
        function testobserved(a, b, ans)
            @test GeneticVariation.DS_DN_NG86(a, b, ncbi_trans_table[1])[1] ≈ ans[1] atol=0.001
            @test GeneticVariation.DS_DN_NG86(a, b, ncbi_trans_table[1])[2] ≈ ans[2] atol=0.001
        end
        answers = [(0.0, 0.0),
                   (0.0, 0.0),
                   (0.0, 1.0),
                   (1.0, 0.0),
                   (1.0, 1.0)]
        for i in 1:length(answers)
            testobserved(codonsA[i], codonsB[i], answers[i])
        end
        testobserved(kmer"TTT", kmer"GAC", (1.0, 2.0))
    end
end
