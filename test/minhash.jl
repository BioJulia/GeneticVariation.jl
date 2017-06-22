@testset "MASH distances" begin
    a = minhash(dna"ATCGCCA-", 4, 3)
    b = minhash(dna"ATCGCCTA", 4, 3)
    @test_approx_eq_eps mash(a, b) 0.2745 1e-3
    @test mash(a, a) == 0
    @test a.sketch == sort(a.sketch)
end
