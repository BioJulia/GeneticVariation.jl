@testset "Proportion distances" begin

    @testset "Biological Sequences" begin
        dnaA = dna"ATCGCCA-"
        dnaB = dna"ATCCCCTA"
        rnaA = rna"ATCGCCA-"
        rnaB = rna"ATCCCCTA"

        @test_approx_eq_eps pdistance(Mismatch, dnaA, dnaB) 0.375 1e-3
        @test_approx_eq_eps pdistance(Mismatch, rnaA, rnaB) 0.375 1e-3

        @test_approx_eq_eps pdistance(Mutated, dnaA, dnaB) 0.2857142857142857 1e-3
        @test_approx_eq_eps pdistance(Mutated, rnaA, rnaB) 0.2857142857142857 1e-3
    end

end
