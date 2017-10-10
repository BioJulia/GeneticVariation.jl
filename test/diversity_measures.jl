@testset "Diversity Measures" begin

    # A simple test to ensure the math works out.
    testP = [0.0 0.006 0.002 0.004;
             0.006 0.0 0.008 0.010;
             0.002 0.008 0.0 0.002;
             0.004 0.010 0.002 0.0;]

    testQ = [0.40, 0.20, 0.20, 0.20]

    @test NL79(testP, testQ) ≈ 0.00352

    # A simple test to ensure the math works out AND the allele
    # frequencies and so on are properly computed from sequence data.
    testSeqs = [dna"AAAACTTTTACCCCCGGGGG",
                dna"AAAACTTTTACCCCCGGGGG",
                dna"AAAACTTTTACCCCCGGGGG",
                dna"AAAACTTTTACCCCCGGGGG",
                dna"AAAAATTTTACCCCCGTGGG",
                dna"AAAAATTTTACCCCCGTGGG",
                dna"AAAACTTTTTCCCCCGTAGG",
                dna"AAAACTTTTTCCCCCGTAGG",
                dna"AAAAATTTTTCCCCCGGAGG",
                dna"AAAAATTTTTCCCCCGGAGG"]

    testP = [0.0 0.1 0.15 0.15;
             0.1 0.0 0.15 0.15;
             0.15 0.15 0.0 0.1;
             0.15 0.15 0.1 0.0;]

    @test NL79(testSeqs) ≈ NL79(testP, testQ)

    testSeqs = [dna"ATAATAAAAAAATAATAAAAAAATAAAAAAAATAAAAAAAA",
                dna"AAAAAAAATAAATAATAAAAAAATAAAAAAAAAAAAAAAAA",
                dna"AAAATAAAAATATAATAAAAAAATATAAAAAAAAAAAAAAA",
                dna"AAAAAAAAAAAATAATAAAAAAATAAATAAATAAAAAAAAA",
                dna"AAAATAAAAAAAATATAAAAAAATAAAAAAAAAAAAAAAAA",
                dna"AAAATAAAAAAAAAATAAAAAAAAAAAAAAAAAAATAAAAA",
                dna"AAAAAATAAAAATAATAAAAAAATAAAAAAAAAAAAAAAAA",
                dna"AAAAAAAAAAAAAAATAAAAAAATAAAAAAAAAAAAAAATA",
                dna"AAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAA",
                dna"AAAAAAAAAAAAAAATAAAAAAATAATAAAAAAAAAAAAAA"]

    @test avg_mut(testSeqs) ≈ 3.888888 atol=10e-5
end
