function apply_prec_ilu0(bilu0, v)
    pv  = zeros( size(v)[1] )
    ccall( (:apply_prec_ILU0, "../extlibs/libmkl.so"), Nothing,
           (Int, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}),
           bilu0.n, bilu0.nzval, bilu0.colptr, bilu0.rowval, v, pv )
    return pv
end
