using LinearAlgebra

function swap_row!(m::Matrix,a,b)
    tmp = m[a,:]
    m[a,:] = m[b,:]
    m[b,:] = tmp
end


function swap_col!(m::Matrix,a,b)
    global column_swaps
    push!(column_swaps,(a,b))
    
    tmp = m[:,a]
    m[:,a] = m[:,b]
    m[:,b] = tmp
end


function nullity(m::Matrix)
    (a,b) = size(m)
    A = m[:,a+1:b]
    return hcat(transpose(-A),I(b-a))
end



function unswap!(m::Matrix)
    global column_swaps

    for x in column_swaps
        tmp = m[:,x[1]]
        m[:,x[1]] = m[:,x[2]]
        m[:,x[2]] = tmp
    end
end



function gaussian_elimination(m::Matrix)
    global column_swaps = Vector{Tuple{Int,Int}}()
    
    m′    = copy(m)
    (a,b) = size(m)
    
    for col in 1:minimum([a,b])
        if !find_pivot_row!(m′,col)
            find_pivot_col!(m′,col)
        end

        if normalize_row!(m′,col)
            eliminate_down!(m′,col)
            eliminate_up!(m′,col)
        end
    end

    return trim_zero_rows(m′)
end


function trim_zero_rows(m::Matrix)
    m′ = copy(m)
    row = size(m′)[1]
    while sum(m′[row,:] .≠ 0)==0
        if row>1
            row -= 1
            m′ = m′[1:row,:]
        else
            break
        end
    end
    return m′
end


function find_pivot_col!(m::Matrix,col)
    col′ = col
    while m[col,col] == 0
        col′ += 1
        if col′ > size(m)[2]
            return false
        end

        if m[col,col′] ≠ 0
            swap_col!(m,col,col′)
        end
    end
    return true
end



function find_pivot_row!(m::Matrix,col)
    row = col + 1
    while m[col,col] == 0
        if row > size(m)[1]
            break
        else
            if m[row,col] ≠ 0
                swap_row!(m,row,col)
            else
                row += 1
            end
        end
    end
    return m[col,col] ≠ 0
end


function eliminate_down!(m::Matrix, col)
    row = col+1
    
    while row <= size(m)[1]
        if m[row,col] ≠ 0
            m[row,:] -= m[row,col].*m[col,:]
        end
        row += 1        
    end
end


function eliminate_up!(m::Matrix, col)
    row = col-1
    
    while row > 0
        if m[row,col] ≠ 0
            m[row,:] -= m[row,col].*m[col,:]
        end
        row -= 1
    end        
end


function normalize_row!(m,col)
    if m[col,col] ≠ 0
        m[col,:] *= 1/m[col,col]
        return true    
    else
        return false
    end
end


    


