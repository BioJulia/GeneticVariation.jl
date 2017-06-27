@testset "NG86" begin
    codons = [kmer"ATG",
              kmer"AAA",
              kmer"CCC",
              kmer"GGG",
              kmer"TTT"]

    @testset "expected" begin
        n_ans = [3.0, 2.666, 2.0, 2.0, 2.666]
        s_ans = [0.0, 0.333, 1.0, 1.0, 0.333]
        for i in 1:endof(codons)
            cdn = codons[i]
            @test_approx_eq_eps GeneticVariation.expected(NG86, cdn, 1.0, ncbi_trans_table[1])[1] s_ans[i] 1e-3
            @test_approx_eq_eps GeneticVariation.expected(NG86, cdn, 1.0, ncbi_trans_table[1])[2] n_ans[i] 1e-3
        end
    end

    @testset "observed" begin
        function testobserved(a, b, ans)
            @test_approx_eq_eps GeneticVariation.observed(NG86, a, b, ncbi_trans_table[1])[1] ans[1] 1e-3
            @test_approx_eq_eps GeneticVariation.observed(NG86, a, b, ncbi_trans_table[1])[2] ans[2] 1e-3
        end
        testobserved(kmer"CCC", kmer"CGC", (0.0, 1.0))
        testobserved(kmer"GGG", kmer"GGC", (1.0, 0.0))
        testobserved(kmer"TTT", kmer"TAC", (1.0, 1.0))
        testobserved(kmer"TTT", kmer"GAC", (1.0, 2.0))
    end
end
