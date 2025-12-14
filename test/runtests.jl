using BlockCodes
using Test

@testset "BlockCodes.jl" begin
    @testset begin
        f=GF(4)
        v=rand(f,3)
        H=rand(f,3,3)
        size(f)
        code=BlockCode(f,H)
        u = sourceword(code)
        c = encode(u,code)
    end
    @testset begin
        println("Verifying simple hamming code for single error correction.")
        Q = rand([2,3,4,5,7,8])
        println("using GF($(Q))")
        code=hamming_code(GF(Q),3)
        u=sourceword(code)
        c=encode(u,code)
        table=syndrome_table(code)
        nonzero_symbols = setdiff(code.field.elements,code.field.F[0])
        for i in 1:length(c)
            for s in nonzero_symbols
                r=copy(c)
                r[i] += s
                pos = table[syndrome(r,code)]
                if pos[1] == i
                    r[pos[1]] += pos[2]
                    if sum(syndrome(r,code))==0
                        println("error corrected at $(pos)")
                    else
                        println("failed correction at $(pos)")
                    end
                else
                    println("wrong error pos $(pos), should be $(i),$(s)")
                end
            end
        end
        @testset begin
            println("Verifying BER for simple binary hamming code.")
            field = GF(2)
            code = hamming_code(field,3)
            m = PSK(2)
            chan = ComplexAWGN(5.0,Q,code.K/code.N,m.Es)
            ec,uc = simulate(code,m,chan,SyndromeDecoder(code))
            pe = 0.0286
            pc = 1-((1-pe)^7 + 7*pe*(1-pe)^6)
            println("Uncoded: expected $(pe), got BER=$(uc.BER)")
            println("Coded:   expected $(pc), got BER=$(ec.BER)")
        end
    end
end
