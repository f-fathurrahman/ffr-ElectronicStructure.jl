function Kprec( pw_grid::PWGrid, psi::Array{ComplexF64,2} )
    Ncol = size(psi)[2]
    G2 = pw_grid.G2
    Npoints = pw_grid.Npoints
    Kpsi = zeros( ComplexF64, size(psi) )
    for ic = 1:Ncol
        for ip = 1:Npoints
          Kpsi[ip,ic] = psi[ip,ic] / ( 1.0 + G2[ip] )
        end
    end
    return Kpsi
end

function Kprec(psi)
  return psi
end
