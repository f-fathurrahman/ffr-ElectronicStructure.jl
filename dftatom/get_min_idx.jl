function get_min_idx( Nmax::Int64, y )
    k = Nmax
    while abs(y[k-1]) < abs(y[k])
        k =  k - 1
        if k == 1 # last index
            break
        end
    end
    k = k - 1
    return k
end