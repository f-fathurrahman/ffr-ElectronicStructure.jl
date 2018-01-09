function norm_matrix_induced(A::Array{Float64,2})
    N = size(A)[1]
    # FIXME no check for matrix form

    # unit-norm vector
    d = 1/sqrt(N)
    v1 = ones(N)*d
    #
    v = A*v1
    return norm(v)
end