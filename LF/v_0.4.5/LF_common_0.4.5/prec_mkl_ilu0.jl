function prec_mkl_ilu0( A::SparseMatrixCSC{Float64,Int64} )
  NNZ = size(A.nzval)[1]
  #println("NNZ = ", NNZ)
  bilu0 = zeros(NNZ)
  #for i = 1:size(A.colptr)[1]
  #  @printf("%8d %8d\n", i, A.colptr[i] )
  #end
  #for i = 1:size(A.rowval)[1]
  #  @printf("%8d %8d %18.10f\n", i, A.rowval[i], A.nzval[i] )
  #end
  ccall( (:mkl_ilu0, "../extlibs/libmkl.so"), Void,
         (Int, Int, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Float64}),
         A.n, NNZ, A.nzval, A.colptr, A.rowval, bilu0 )
  return SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, bilu0)
end
