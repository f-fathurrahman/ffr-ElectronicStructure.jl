"""
Test orthonormalization of wavefunction in real space
"""
function test_ortho( dVol, Nstates, psi )
  print("\nTest norm:\n")
  print( "===========================\n")
  for is = 1:Nstates
   print("State: #$(is): $(dot( psi[:,is], psi[:,is] ) * dVol)\n")
  end
  print("\nTest ortho w.r.t state #1:\n")
  print( "===========================\n")
  for is = 2:Nstates
   print("State: #$(is): $(dot( psi[:,is], psi[:,1] ) * dVol)\n")
  end
  print("\n")
end
