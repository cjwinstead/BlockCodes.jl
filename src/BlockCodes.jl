module BlockCodes

using Parameters, Primes, GaloisFields, LinearAlgebra, DataFrames, CSV


include("GF.jl")
include("gaussian_elimination.jl")

include("Modulations.jl")
include("ChannelModels.jl")
include("ErrorCounter.jl")

"""
   BlockCode struct 

Represent block error-correcting code over GF(Q), 
where Q is a power of 2.

BlockCode(
    field::GF
    H::Matrix
    G::Matrix
    N::Int
    K::Int
    R::Float64
    sym_node_neighbors::Vector{Vector}
    check_node_neighbors::Vector{Vector}
    M::Int
)
"""
@with_kw struct BlockCode
    field::GF
    H::Matrix
    G::Matrix
    N::Int
    K::Int
    R::Float64
    sym_node_neighbors::Vector{Vector}
    check_node_neighbors::Vector{Vector}
    M::Int
end


function BlockCode(field::GF,filename::String)
    m = Matrix(DataFrame(CSV.File(filename;header=0)))
    H = map.(x -> field.elements[x+1], m)
    return BlockCode(field,H)
end



function BlockCode(field::GF,H::Matrix)
    
    M,N = size(H)
    
    check_node_neighbors = Vector{Vector{Int}}()
    sym_node_neighbors   = Vector{Vector{Int}}()

    H′ = gaussian_elimination(H)
    G  = nullity(H′)
    if minimum(size(G))==0
        G = H′
    else
        unswap!(G)
    end
    
    K = size(G)[1]
    R = K/N
    
    for i in 1:M
        this_checks_neighbors = Vector{Int}()
        for j in 1:N
            if H[i,j] ≠ 0
                push!(this_checks_neighbors,j)
            end
        end
        push!(check_node_neighbors,this_checks_neighbors)
    end

    for j in 1:N
        this_syms_neighbors = Vector{Int}()
        for i in 1:M
            if H[i,j] ≠ 0
                push!(this_syms_neighbors,i)
            end
        end
        push!(sym_node_neighbors,this_syms_neighbors)
    end
    
    return BlockCode(
        field,
        H,
        G,
        N,
        K,
        R,
        sym_node_neighbors,
        check_node_neighbors,
        M
    )

end


"""
   syndrome(v::Vector{GF(Q)}, H::Matrix{GF(Q)})

Compute syndrome vector (i.e. parity-checks) for symbols v
with parity-check matrix H.

Returns Vector{GF(Q)} of length M
"""
function syndrome(v::Vector,H::Matrix)
    return H*v
end

"""
   syndrome(v::Vector{GF(Q)}, code::BlockCode)

Compute syndrome vector (i.e. parity-checks) for symbols v
with parity-check matrix H defined by `code`.

Returns Vector{GF(Q)} of length M
"""
function syndrome(v::Vector,code::BlockCode)
    return code.H*v
end


"""
    encode(u::Vector{GF(Q)},code::BlockCode)

Encode source word u with generator matrix defined by `code`.

Returns Vector{GF(Q)} of length N
"""
function encode(u::Vector,code::BlockCode)
    return transpose(code.G)*u
end


"""
    sourceword(code::BlockCode)

Return Vector{GF(Q)} of random symbols with length K
"""
function sourceword(code::BlockCode)
    return rand(code,code.K)
end


function rand(code::BlockCode,a::Int64)
    return rand(code.field,a)
end

function transmit(v::Vector,modulation::Modulation,chan::ChannelModel)
    return modulate(v,modulation) .+ noise_samples(chan,length(v))
end

function transmit(v::Vector,chan::ChannelModel)
    return v .+ noise_samples(chan,length(v))
end


function isuniquecolumn(f::GF,v,H)
    if size(H)[1] == 0
        return true
    end
    
    for i in 1:size(H)[2]
        for e in f
            if (e ≠ 0)
#                println("$(i): $(e) .* $(v) vs $(H[:,i])")
 #               println("      eq? $(e.*v .== H[:,i])")
                if prod(e.*v .== H[:,i])
                    return false
                end
            end
        end
    end
    return true
end

function hamming_code(f::GF,M::Int64)
    N = f.Q^M - 1
    H = f.F[]
   
    for i in 1:N        
        col = to_gf(digits(i;base=f.Q,pad=M),f)
        if (i==1)
            H = hcat(Vector{f.F}(col))
        elseif isuniquecolumn(f,col,H)
            H = hcat(H,Vector{f.F}(col))
        end
    end
    return BlockCode(f,H)
end


function syndrome_table(code::BlockCode)
    f = code.field
    table = Dict{Vector{f.F},Tuple{Int64,f.F}}()
    for i in 1:code.N
        for symbol in code.field
            if symbol ≠ 0
                c = zeros(f.F,code.N)
                c[i] = symbol
                S = syndrome(c,code)
                if sum(S .≠ 0)>0
                    if haskey(table,S)
                        println("i=$(i) sym=$(symbol) key $(S) already assigned to $(table[S])")
                    else
                        table[S] = (i,-symbol)
                    end
                end
            end
        end
    end

    return table
end


function ComplexAWGN(Eb_N₀::Float64,code::BlockCode,modulation::Modulation)
    R = code.R# * log2(code.field.Q)
    Q = code.field.Q
    Es = modulation.Es
    
    return ComplexAWGN(Eb_N₀,Q,R,Es)
end

include("Decoder.jl")


function transaction(code::BlockCode,m::Modulation,chan::ChannelModel)
    u = sourceword(code)
    c = encode(u,code)
    x = modulate(c,m)
    y = transmit(x,chan)
    r = hard_decision(y,m)

    return u,c,x,y,r
end


function simulate(code::BlockCode,
                  m::Modulation,
                  chan::ChannelModel,
                  dec::Decoder;
                  maxwords=1e6,
                  errwords=150)
    
    ec = ErrorCounter(code.field)
    uc = ErrorCounter(code.field)
    while (ec.word_errors < errwords) & (ec.total_words < maxwords)
        u,c,x,y,r = transaction(code,m,chan)
        d = decode(y,r,dec)
        counterrors!(ec,c,d)
        counterrors!(uc,c,r)
    end
    return ec,uc
end

         

export GF
export rand
export size
export iterate
export getindex
export to_binary
export BlockCode
export sourceword
export encode
export syndrome
export hamming_code
export syndrome_table
export gaussian_elimination
export nullity
export swap_row!
export swap_col!
export unswap!
export trim_zero_rows
export find_pivot_col!
export find_pivot_row!
export eliminate_down!
export eliminate_up!
export normalize_row!
export Modulation
export RealModulation
export ComplexModulation
export ASK
export PSK
export QAM
export modulate
export hard_decision
export transmit
export ChannelModel
export RealChannel
export ComplexChannel
export ComplexAWGN
export RealAWGN
export AWGN
export noise_samples
export transaction
export ErrorCounter
export counterrors!
export simulate
export Decoder
export SyndromeDecoder
export isuniquecolumn
export to_gf

end
