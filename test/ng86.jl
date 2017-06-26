@testset "NG86" begin
    codons = [kmer"ATG",
              kmer"AAA",
              kmer"CCC",
              kmer"GGG",
              kmer"TTT",
              kmer"TAA"]

    @testset "expected" begin
        n_ans = [3.0, 3.666, 2.0, 2.0, 2.666, 2.333]
        s_ans = [0.0, 0.333, 1.0, 1.0, 0.333, 0.666]
        for i in 1:endof(codons)
            cdn = codons[i]
            @test_approx_eq_eps GeneticVariation.expected(cdn)[1] s_ans[i] 1e-3
            @test_approx_eq_eps GeneticVariation.expected(cdn)[2] n_ans[i] 1e-3
        end
    end

end
