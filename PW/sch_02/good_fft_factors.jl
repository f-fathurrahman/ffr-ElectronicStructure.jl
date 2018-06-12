function is_good_factors( N::Int )

    good_factors = [2, 3, 5, 7]

    for k in keys( factor(N) )
        is_in = k in good_factors
        if !is_in
            return false
        end
    end
    return true

end
