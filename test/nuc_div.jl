@testset "Nucleotide Diversity" begin

    # A simple test to ensure the math works out.
    testP = [0.0 0.006 0.002 0.004;
             0.006 0.0 0.008 0.010;
             0.002 0.008 0.0 0.002;
             0.004 0.010 0.002 0.0;]

    testQ = [0.40, 0.20, 0.20, 0.20]

    @test nuc_div(testP, testQ) ≈ 0.00352


end
