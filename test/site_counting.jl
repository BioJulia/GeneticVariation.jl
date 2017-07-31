@testset "Site counting" begin

    dna_alphabets = (DNAAlphabet{4}, DNAAlphabet{2})
    rna_alphabets = (RNAAlphabet{4}, RNAAlphabet{2})

    @testset "Specific count methods" begin

        function generate_possibilities_tester(::Type{A}) where A <: NucAlphs
            symbols = alphabet(A)
            arra = Vector{eltype(A)}()
            arrb = Vector{eltype(A)}()
            for i in 1:length(symbols), j in i:length(symbols)
                push!(arra, symbols[i])
                push!(arrb, symbols[j])
            end
            return BioSequence{A}(arra), BioSequence{A}(arrb)
        end

        @testset "2 bit" begin
            dnaA, dnaB = generate_possibilities_tester(DNAAlphabet{2})
            rnaA, rnaB = generate_possibilities_tester(RNAAlphabet{2})

            @test count(Conserved, dnaA, dnaB) == count(Conserved, dnaB, dnaA) == 4
            @test count(Mutated, dnaA, dnaB) == count(Mutated, dnaB, dnaA) == 6

            @test count(Conserved, rnaA, rnaB) == count(Conserved, rnaB, rnaA) == 4
            @test count(Mutated, rnaA, rnaB) == count(Mutated, rnaB, rnaA) == 6
        end

        @testset "4 bit" begin
            dnaA, dnaB = generate_possibilities_tester(DNAAlphabet{4})
            rnaA, rnaB = generate_possibilities_tester(RNAAlphabet{4})

            @test count(Conserved, dnaA, dnaB) == count(Conserved, dnaB, dnaA) == (4, 10)
            @test count(Mutated, dnaA, dnaB) == count(Mutated, dnaB, dnaA) == (6, 10)

            @test count(Conserved, rnaA, rnaB) == count(Conserved, rnaB, rnaA) == (4, 10)
            @test count(Mutated, rnaA, rnaB) == count(Mutated, rnaB, rnaA) == (6, 10)
        end
    end

    @testset "Randomized tests" begin
        # A test counting function which is naive.
        @inline function issite(::Type{Conserved}, a::BioSequence, b::BioSequence, idx)
            return a[idx] == b[idx]
        end
        @inline function issite(::Type{Mutated}, a::BioSequence, b::BioSequence, idx)
            return a[idx] != b[idx]
        end
        @inline function testcount2(::Type{P}, a::BioSequence, b::BioSequence) where P <: BioSequences.Position
            k = 0
            @inbounds for idx in 1:min(endof(a), endof(b))
                k += issite(P, a, b, idx)
            end
            return k
        end
        @inline function testcount(::Type{P}, a::BioSequence, b::BioSequence) where P <: BioSequences.Position
            k, c = 0, 0
            @inbounds for idx in 1:min(endof(a), endof(b))
                isvalid = !(isgap(a[idx]) || isgap(b[idx])) && !(isambiguous(a[idx]) || isambiguous(b[idx]))
                k += issite(P, a, b, idx) && isvalid
                c += isvalid
            end
            return k, c
        end
        function testcounting(::Type{S}, a::BioSequence, b::BioSequence) where S <: Site
            @test count(S, a, b) == count(S, b, a) == testcount(S, a, b)
        end
        function testcounting2(::Type{S}, a::BioSequence, b::BioSequence) where S <: Site
            @test count(S, a, b) == count(S, b, a) == testcount2(S, a, b)
        end
        function testforencs(a::Int, b::Int, subset::Bool)
            for alphabet in (DNAAlphabet, RNAAlphabet)
                for _ in  1:10
                    seqA = random_seq(alphabet{a}, rand(10:100))
                    seqB = random_seq(alphabet{b}, rand(10:100))
                    sa = seqA
                    sb = seqB
                    if subset
                        intA = random_interval(1, length(seqA))
                        intB = random_interval(1, length(seqB))
                        subA = seqA[intA]
                        subB = seqB[intB]
                        sa = subA
                        sb = subB
                    end
                    if a == 2 && b == 2
                        testcounting2(Conserved, sa, sb)
                        testcounting2(Mutated, sa, sb)
                    else
                        testcounting(Conserved, sa, sb)
                        testcounting(Mutated, sa, sb)
                    end
                end
            end
        end

        @testset "Testing 4-bit seq pairs" begin
            @testset "Full random sequences" begin
                testforencs(4, 4, false)
            end
            @testset "Subset random sequences" begin
                testforencs(4, 4, true)
            end
        end

        @testset "Testing 2-bit seq pairs" begin
            @testset "Full random sequences" begin
                testforencs(2, 2, false)
            end
            @testset "Subset random sequences" begin
                testforencs(2, 2, true)
            end
        end

        @testset "Testing mixed bit seq pairs" begin
            @testset "Full random sequences" begin
                testforencs(4, 2, false)
                testforencs(2, 4, false)
            end
            @testset "Subset random sequences" begin
                testforencs(4, 2, true)
                testforencs(2, 4, true)
            end
        end
    end

    @testset "Pairwise methods" begin
        @testset "4-bit encoded sequences" begin
            dnas = [dna"ATCGCCA-", dna"ATCGCCTA", dna"ATCGCCT-", dna"GTCGCCTA"]
            rnas = [rna"AUCGCCA-", rna"AUCGCCUA", rna"AUCGCCU-", rna"GUCGCCUA"]

            answer_mutated = [(0,0) (1,7) (1,7) (2,7);
                              (1,7) (0,0) (0,7) (1,8);
                              (1,7) (0,7) (0,0) (1,7);
                              (2,7) (1,8) (1,7) (0,0)]

            answer_conserved = [(0,0) (6,7) (6,7) (5,7);
                                (6,7) (0,0) (7,7) (7,8);
                                (6,7) (7,7) (0,0) (6,7);
                                (5,7) (7,8) (6,7) (0,0)]
            for i in (dnas, rnas)
                @test count_pairwise(Mutated, i...) == answer_mutated
                @test count_pairwise(Conserved, i...) == answer_conserved
            end
        end
        @testset "2-bit encoded sequences" begin
            dnas = [BioSequence{DNAAlphabet{2}}("ATCGCCAC"),
                    BioSequence{DNAAlphabet{2}}("ATCGCCTA"),
                    BioSequence{DNAAlphabet{2}}("ATCGCCTT"),
                    BioSequence{DNAAlphabet{2}}("GTCGCCTA")]
            rnas = [BioSequence{RNAAlphabet{2}}("AUCGCCAC"),
                    BioSequence{RNAAlphabet{2}}("AUCGCCUA"),
                    BioSequence{RNAAlphabet{2}}("AUCGCCUU"),
                    BioSequence{RNAAlphabet{2}}("GUCGCCUA")]
            answer_mutated = [0 2 2 3;
                              2 0 1 1;
                              2 1 0 2;
                              3 1 2 0]
            answer_conserved = [0 6 6 5;
                                6 0 7 7;
                                6 7 0 6;
                                5 7 6 0]
            for i in (dnas, rnas)
                @test count_pairwise(Mutated, i...) == answer_mutated
                @test count_pairwise(Conserved, i...) == answer_conserved
            end
        end
    end

    @testset "Windowed methods" begin
        @testset "4-bit encoded sequences" begin
            dnaA = dna"ATCGCCA-M"
            dnaB = dna"ATCGCCTAA"
            rnaA = rna"AUCGCCA-M"
            rnaB = rna"AUCGCCUAA"

            for seqs in ((dnaA, dnaB), (rnaA, rnaB))
                @test count(Conserved, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, (3,3)),
                                                                   IntervalValue(2, 4, (3,3)),
                                                                   IntervalValue(3, 5, (3,3)),
                                                                   IntervalValue(4, 6, (3,3)),
                                                                   IntervalValue(5, 7, (2,3)),
                                                                   IntervalValue(6, 8, (1,2)),
                                                                   IntervalValue(7, 9, (0,1))]
                @test count(Mutated, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, (0,3)),
                                                                 IntervalValue(2, 4, (0,3)),
                                                                 IntervalValue(3, 5, (0,3)),
                                                                 IntervalValue(4, 6, (0,3)),
                                                                 IntervalValue(5, 7, (1,3)),
                                                                 IntervalValue(6, 8, (1,2)),
                                                                 IntervalValue(7, 9, (1,1))]
            end
        end

        @testset "2-bit encoded sequences" begin
            dnaA = BioSequence{DNAAlphabet{2}}("ATCGCCATT")
            dnaB = BioSequence{DNAAlphabet{2}}("ATCGCCTAA")
            rnaA = BioSequence{RNAAlphabet{2}}("AUCGCCAUU")
            rnaB = BioSequence{RNAAlphabet{2}}("AUCGCCUAA")

            for seqs in ((dnaA, dnaB), (rnaA, rnaB))
                @test count(Conserved, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 3),
                                                                   IntervalValue(2, 4, 3),
                                                                   IntervalValue(3, 5, 3),
                                                                   IntervalValue(4, 6, 3),
                                                                   IntervalValue(5, 7, 2),
                                                                   IntervalValue(6, 8, 1),
                                                                   IntervalValue(7, 9, 0)]
                @test count(Mutated, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, 0),
                                                                 IntervalValue(2, 4, 0),
                                                                 IntervalValue(3, 5, 0),
                                                                 IntervalValue(4, 6, 0),
                                                                 IntervalValue(5, 7, 1),
                                                                 IntervalValue(6, 8, 2),
                                                                 IntervalValue(7, 9, 3)]
            end
        end

        @testset "Mixed encodings" begin
            dnaA = dna"ATCGCCA-M"
            dnaB = BioSequence{DNAAlphabet{2}}("ATCGCCTAA")
            rnaA = rna"AUCGCCA-M"
            rnaB = BioSequence{RNAAlphabet{2}}("AUCGCCUAA")

            for seqs in ((dnaA, dnaB), (rnaA, rnaB))
                @test count(Conserved, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, (3,3)),
                                                                   IntervalValue(2, 4, (3,3)),
                                                                   IntervalValue(3, 5, (3,3)),
                                                                   IntervalValue(4, 6, (3,3)),
                                                                   IntervalValue(5, 7, (2,3)),
                                                                   IntervalValue(6, 8, (1,2)),
                                                                   IntervalValue(7, 9, (0,1))]
                @test count(Mutated, seqs[1], seqs[2], 3, 1) == [IntervalValue(1, 3, (0,3)),
                                                                 IntervalValue(2, 4, (0,3)),
                                                                 IntervalValue(3, 5, (0,3)),
                                                                 IntervalValue(4, 6, (0,3)),
                                                                 IntervalValue(5, 7, (1,3)),
                                                                 IntervalValue(6, 8, (1,2)),
                                                                 IntervalValue(7, 9, (1,1))]
            end
        end
    end
end
