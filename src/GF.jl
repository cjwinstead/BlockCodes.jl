import Base.rand
import Base.size

𝔽₂ = @GaloisField 2

@with_kw struct GF
    Q::Int64 
    F
    β
    elements::Vector
    order::Dict
end

"""
    GF(q::Int64)

Initialize GF (Galois Field) struct with field order Q.
"""
function GF(q::Int64)
    
    if isprime(q)
        F = @GaloisField q^1
        β = nothing
    elseif length(factor(Set,q))==1

        F = @GaloisField! q β
    
    else
        return nothing
    end
    
    order,elements = define_symbol_order(F)
    return GF(Q=q,
              β=β,
              F=F,
              elements=elements,
              order=order
              )
end

"""
    rand(f::GF,a::Int64,b::Int64)

Return randomly selected elements from GF(Q).
"""
function rand(f::GF,a::Int64,b::Int64)
    return rand(f.F,a,b)
end


"""
    rand(f::GF,a::Int64)

Return randomly selected elements from GF(Q).
"""
function rand(f::GF,a::Int64)
    return rand(f.F,a)
end



"""
    size(f::GF)

Returns number of elements in field f.
"""
function size(f::GF)
    return f.Q
end


"""
    iterate(f::GF)

Iterate over field elements in f.

Returns element,state
"""
function Base.iterate(f::GF)

    if f.Q == 0
        return nothing,nothing
    else
        state = Int(1)
        return f.elements[1],state+1
    end
end


"""
    iterate(f::GF,state)

Iterate over field elements in f.

Returns element,state
"""
function Base.iterate(f::GF,state)
    if isnothing(state)
        return nothing
    elseif f.Q == state
        return f.elements[state],nothing
    else
        return f.elements[state],state+1
    end
end


"""
    getindex(f::GF,index::Int)

Return element from f at specified index.
"""
function Base.getindex(f::GF,index::Int64)
    if index > 0 && index <= f.Q
        return f.elements[index]
    else
        return nothing
    end
end



function define_symbol_order(F)
    
    order = Dict{F,F}()
    elements = Vector{F}()
    
    let sym = 0
        for x in F
            push!(elements,x)
            if (x ≠ 0)            
                order[sym] = x
                sym = x
            end
        end
    
        order[sym] = 0
    end
    return order,elements
end


"""
    to_binary(v::Vector{GF(Q)},f::GF(Q))

Returns Vector{𝔽₂} representing the elements of v.
"""
function to_binary(v::Vector,field::GF)
    b = Vector{𝔽₂}()
    for x in v
        b=vcat(b,digits(x.n;base=2,pad=ceil(Int,log2(field.Q))))
    end
    return b
end


function to_gf(x::Int64,f::GF)
    return f.elements[x+1]
end

function to_gf(v::Vector{Int64},f::GF)
    return map(x->f.elements[x+1],v)
end
