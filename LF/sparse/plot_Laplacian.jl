const ⊗ = kron

import PyPlot
const plt = PyPlot

function plot_matrix( Nx, Ny, Nz, matrix, postfix )
  plt.clf()
  plt.spy(matrix)
  desc = string(Nx)*"_"*string(Ny)*"_"*string(Nz)
  filename = desc*"_"*postfix*".png"
  plt.savefig(filename, dpi=300)

  println("")
  println(desc*postfix)
  println(sparse(matrix))
end

function test_main( Nx, Ny, Nz )

  D2x = rand(Nx,Nx); D2x = 0.5*( D2x + D2x' )
  D2y = rand(Ny,Ny); D2y = 0.5*( D2y + D2y' )
  D2z = rand(Nz,Nz); D2z = 0.5*( D2z + D2z' )

  term1 = D2x ⊗ eye(Ny,Ny) ⊗ eye(Nz,Nz)
  term2 = eye(Nx,Nx) ⊗ D2y ⊗ eye(Nz,Nz)
  term3 = eye(Nx,Nx) ⊗ eye(Ny,Ny) ⊗ D2z
  Laplacian = term1 + term2 + term3

  plot_matrix( Nx, Ny, Nz, Laplacian, "Laplacian")
  plot_matrix( Nx, Ny, Nz, term1, "term1")
  plot_matrix( Nx, Ny, Nz, term2, "term2")
  plot_matrix( Nx, Ny, Nz, term3, "term3")

end

#test_main( 5, 3, 3 )
#test_main( 3, 5, 3 )
#test_main( 3, 3, 5 )

test_main( 3, 3, 3 )
