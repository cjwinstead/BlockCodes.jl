

abstract type ChannelModel   end
abstract type ComplexChannel <: ChannelModel end
abstract type RealChannel    <: ChannelModel end

@with_kw struct ComplexAWGN <: ComplexChannel
    Eb_N₀::Float64
    Q::Int64
    R::Float64
    Es::Float64
    N₀::Float64
    σ₀::Float64
end

function ComplexAWGN(Eb_N₀::Float64,Q::Int64,R::Float64,Es::Float64)
    Eb_N₀_linear = 10 .^(Eb_N₀/10)
    N₀ = Es/(log2(Q)*R*Eb_N₀_linear) 
    σ₀ = sqrt(N₀/2)

    return ComplexAWGN(Eb_N₀,Q,R,Es,N₀,σ₀)
end




function noise_samples(chan::ComplexAWGN,N::Int)
    chan.σ₀.*(randn(N) .+ 1im.*randn(N))
end



