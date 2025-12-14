
@with_kw mutable struct ErrorCounter
    field::GF
    total_words::Integer = 0
    total_symbols::Integer = 0
    total_bits::Integer = 0
    word_errors::Integer = 0
    symbol_errors::Integer = 0
    bit_errors::Integer = 0
    SER::Float64 = 0
    BER::Float64 = 0
    WER::Float64 = 0
end

function ErrorCounter(f::GF)
    return ErrorCounter(f,0,0,0,
                        0,0,0,
                        0,0,0)
end


function counterrors!(cntr::ErrorCounter,r::Vector,c::Vector)
    @unpack_ErrorCounter cntr
    
    new_symbol_errors = sum(r .≠ c)
    if (new_symbol_errors > 0)
        symbol_errors  += new_symbol_errors
        new_bit_errors  = sum(to_binary(r,field) .≠ to_binary(c,field))
        bit_errors     += new_bit_errors
        word_errors    += 1
    end
    total_words   += 1
    total_symbols += length(r)
    total_bits    += length(r) * log2(field.Q)
    SER = symbol_errors/total_symbols
    BER = bit_errors/total_bits
    WER = word_errors/total_words

    @pack_ErrorCounter! cntr
end
