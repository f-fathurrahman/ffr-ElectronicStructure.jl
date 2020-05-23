using Printf
using SparseArrays

import PyPlot
const plt = PyPlot

include("../LF_common/m_LF1d.jl")
include("../LF_common/m_LF3d.jl")

function test_main( NN::Array{Int64} )
  AA = [0.0, 0.0, 0.0]
  BB = [16.0, 16.0, 16.0]
  LF = init_LF3d_c( NN, AA, BB )

  Npoints = prod(NN)
  println("Step 1: ")
  
  NNZ = ccall( (:get_nnz_laplacian3d_, "./Laplace3d_sparse.so"), Int64,
         (Ptr{Int64}, Ref{Int64}, Ptr{Int64},
          Ptr{Float64}, Ptr{Float64}, Ptr{Float64} ),
         NN, Npoints, LF.lin2xyz, LF.LFx.D2jl, LF.LFy.D2jl, LF.LFz.D2jl )

  #
  row    = zeros(Int64,NNZ)
  column = zeros(Int64,NNZ)
  values = zeros(Float64,NNZ)
  #
  println("Step 2: ")

  ccall( (:init_laplacian3d_, "./Laplace3d_sparse.so"), Nothing,
         ( Ref{Int64}, Ptr{Int64}, Ref{Int64}, Ptr{Int64},
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Int64}, Ptr{Int64}, Ptr{Float64} ),
         NNZ, NN, Npoints, LF.lin2xyz, LF.LFx.D2jl, LF.LFy.D2jl, LF.LFz.D2jl,
         row, column, values )

  #
  #for ii = 1:NNZ
  #  @printf("%5d %5d %18.10f\n", row[ii], column[ii], values[ii])
  #end
  println("Step 3: ")

  Laplacian = sparse( row, column, values )

  plt.clf()
  plt.spy(Matrix(Laplacian))
  filename = "Laplacian_"*string(NN[1])*"_"*string(NN[2])*"_"*string(NN[3])*".pdf"
  plt.savefig(filename)

end

test_main([3,5,5])
#test_main([20,20,20])
#test_main([30,30,30])
#test_main([40,40,40])
#test_main([50,50,50])
#for ii = 10:10:60
#  @time test_main( ii*ones(Int64,3) )
#end
