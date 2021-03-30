function get_n_nodes(N, y)
    Nnodes = 0
    last_sign = Int( sign(y[1]) )
    last_i = -1
    for i = 2:N
        isy = Int( sign(y[i]) )
        if isy == -last_sign
            last_sign = isy
            last_i = i - 1
            Nnodes = Nnodes + 1
        end
    end
    return Nnodes
end