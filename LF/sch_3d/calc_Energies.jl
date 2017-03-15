function calc_Energies( LF::LF3dGrid, Vpot, psi::Array{Float64,2} )
  #
  Ncol = size(psi)[2]
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
  Etot = 0.0
  Ekin = 0.0
  Epot = 0.0
  #
  for ic = 1:Ncol
    nabla2_v = op_nabla2( LF, psi[:,ic] )
    Ekin = Ekin + -0.5*dot( psi[:,ic], nabla2_v ) * ΔV
    Epot = Epot + sum( Vpot .* psi[:,ic].^2 ) * ΔV
  end
  Etot = Ekin + Epot
  return EnergiesT( Etot, Ekin, Epot )
end


function calc_Energies( LF::LF3dGrid, ∇2::SparseMatrixCSC{Float64,Int64},
                        Vpot, psi::Array{Float64,2} )
  #
  Ncol = size(psi)[2]
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
  Etot = 0.0
  Ekin = 0.0
  Epot = 0.0
  #
  for ic = 1:Ncol
    nabla2_v = ∇2 * psi[:,ic]
    Ekin = Ekin + -0.5*dot( psi[:,ic], nabla2_v ) * ΔV
    Epot = Epot + sum( Vpot .* psi[:,ic].^2 ) * ΔV
  end
  Etot = Ekin + Epot
  return EnergiesT( Etot, Ekin, Epot )
end
