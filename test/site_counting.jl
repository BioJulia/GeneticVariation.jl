@testset "Site counting" begin

    dna_alphabets = (DNAAlphabet{4}, DNAAlphabet{2})
    rna_alphabets = (RNAAlphabet{4}, RNAAlphabet{2})

    @testset "Specific count methods" begin

        function generate_possibilities_tester{A<:NucleicAcidAlphabets}(::Type{A})
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
        @inline function testcount{P<:BioSequences.Position}(::Type{P}, a::BioSequence, b::BioSequence)
            k = 0
            @inbounds for idx in 1:min(endof(a), endof(b))
                k += issite(P, a, b, idx)
            end
            return k
        end

        @inline function issite(::Type{Conserved}, a::BioSequence, b::BioSequence, idx)
            return !(isambiguous(a[idx]) || isambiguous(b[idx])) && (a[idx] == b[idx])
        end
        @inline function issite(::Type{Mutated}, a::BioSequence, b::BioSequence, idx)
            return !(isambiguous(a[idx]) || isambiguous(b[idx])) && (a[idx] != b[idx])
        end

        function testcounting{S<:Site}(::Type{S}, a, b)
            @test count(S, a, b) == count(S, b, a) == testcount(S, a, b)
        end

        function testforencs(a::Int, b::Int, subset::Bool)
            for alphabet in (DNAAlphabet, RNAAlphabet)
                for _ in  1:50
                    #println("TESTING SEQUENCES:")
                    seqA = random_seq(alphabet{a}, rand(10:100))
                    seqB = random_seq(alphabet{b}, rand(10:100))
                    #=
                    println("A seq: ", seqA)
                    println("B seq: ", seqB)
                    println("A length: ", length(seqA))
                    println("B length: ", length(seqB))
                    =#
                    sa = seqA
                    sb = seqB
                    if subset
                        intA = random_interval(1, length(seqA))
                        intB = random_interval(1, length(seqB))
                        subA = seqA[intA]
                        subB = seqB[intB]
                        #=
                        println("A subset: ", intA)
                        println("B subset: ", intB)
                        println("A subseq: ", subA)
                        println("B subseq: ", subB)
                        =#
                        sa = subA
                        sb = subB
                    end
                    testcounting(Conserved, sa, sb)
                    testcounting(Mutated, sa, sb)
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
end
