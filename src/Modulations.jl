
abstract type Modulation end
abstract type RealModulation <: Modulation end
abstract type ComplexModulation <: Modulation end

function modulate(v::Vector,m::Modulation)
    return map(x->m.symbol_map[x],v) 
end


@with_kw struct ASK <: RealModulation
    field
    constellation
    symbol_map
    Es
end

function ASK(Q::Int64)
    field = GF(Q)
    return ASK(field)
end

function ASK(field::GF)    
    symbol_map=Dict{GaloisFields.AbstractGaloisField,Float64}()
    Es=0.0
    mags = 2.0.*(0.5 .+ collect(0:(field.Q/2-1)))
    constellation = (-1).^collect(1:field.Q) .* repeat(mags,inner=2)
    
    for i in 1:field.Q
        symbol_map[field[i]] = constellation[i]
        Es += constellation[i]^2
    end
    return ASK(field,constellation,symbol_map,Es/field.Q)
end


##########################################################

@with_kw struct PSK <: ComplexModulation
    field
    constellation
    symbol_map
    Es
end

function PSK(Q::Int64)
    field = GF(Q)
    return PSK(field)
end

function PSK(field::GF)    
    symbol_map=Dict{GaloisFields.AbstractGaloisField,Complex}()
    Es=1.0

    constellation = round.(Es.*exp.(1im.*pi.*2.0.*collect(0:field.Q-1)./field.Q);digits=2)
    
    for i in 1:field.Q
        symbol_map[field[i]] = constellation[i]
    end
    return PSK(field,constellation,symbol_map,Es)
end

##########################################################

@with_kw struct QAM <: ComplexModulation
    field
    symbol_distance
    constellation
    symbol_map
    Es
end

function QAM(Q::Int64)
    field = GF(Q)
    return QAM(field)
end

function QAM(field::GF)
    return QAM(field,1.0)
end

function QAM(field::GF, symbol_distance::Float64)
             
    T = sqrt(field.Q)
    i = 1

    constellation = Vector{Complex}()
    symbol_map    = Dict{field.F}{Complex}()

    
    for j in collect((-T/2+0.5):(T/2))
        for k in collect((-T/2+0.5):(T/2))
            push!(constellation,symbol_distance*((j) + (k)*1im))
            symbol_map[field.elements[i]] = constellation[i]
            i += 1
        end
    end

    Es = real(constellation'*constellation/field.Q)

    return QAM(
        field,
        symbol_distance,
        constellation,
        symbol_map,
        Es
    )

end



function hard_decision(v::Vector,m::Modulation)
    d = Vector{m.field.F}(undef,length(v))
                     
    for i in 1:length(v)
        δ        = abs.(v[i] .- m.constellation)
        δmin,idx = findmin(δ)
        d[i]     = findfirst(y->y==m.constellation[idx],m.symbol_map)
    end
    return d
end


