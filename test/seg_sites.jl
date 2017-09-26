@testset "Segregating sites" begin

    S = [dna"aatcga", dna"aaa", dna"tatcg", dna"catcgac", dna"aatc"]

    @test nsegregating(S)[1] == 2 && nsegregating(S)[2] == 3
end
