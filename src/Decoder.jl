

abstract type Decoder end

@with_kw struct SyndromeDecoder <: Decoder
    syndrome_table
    code
end

function SyndromeDecoder(code::BlockCode)
    return SyndromeDecoder(syndrome_table(code),code)
end

function decode(y::Vector,r::Vector,dec::SyndromeDecoder)
    return decode(r,dec)
end


function decode(r::Vector,dec::SyndromeDecoder)
    s = syndrome(r,dec.code)
    if sum(s) == 0
        return r
    else
        d = copy(r)
        (i,e) = dec.syndrome_table[s]
        d[i] += e
        
        return d
    end
end
