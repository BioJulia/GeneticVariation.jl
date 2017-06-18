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

        for alphset in (dna_alphabets, rna_alphabets)
            # Answers to these tests were worked out manually to verify
            # counting is working correctly.
            # seqA and seqB contain all possible observations of sites.

            # 4 bit encoded sequences
            seqA, seqB = generate_possibilities_tester(alphset[1])
            @test count(Conserved, seqA, seqB) == count(Conserved, seqB, seqA) == length(alphabet(alphset[1]))
            @test count(Mutated, seqA, seqB) == count(Mutated, seqB, seqA) == (length(seqA) - length(alphabet(alphset[1])))

            # 2 bit encoded sequences
            seqA, seqB = generate_possibilities_tester(alphset[2])
            @test count(Conserved, seqA, seqB) == count(Conserved, seqB, seqA) == length(alphabet(alphset[2]))
            @test count(Mutated, seqA, seqB) == count(Mutated, seqB, seqA) == (length(seqA) - length(alphabet(alphset[2])))
        end
    end
end
