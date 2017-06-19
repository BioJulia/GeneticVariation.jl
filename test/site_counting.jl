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

            @test count(Conserved, dnaA, dnaB) == count(Conserved, dnaB, dnaA) == 0
            @test count(Mutated, dnaA, dnaB) == count(Mutated, dnaB, dnaA) == 0

            @test count(Conserved, rnaA, rnaB) == count(Conserved, rnaB, rnaA) == 0
            @test count(Mutated, rnaA, rnaB) == count(Mutated, rnaB, rnaA) == 0
        end
    end
end
