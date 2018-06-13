function Kprec( pw_grid::PWGrid, psi::Array{ComplexF64,2} )
    Nstates = size(psi)[2]
    G2 = pw_grid.G2
    Npoints = pw_grid.Npoints
    Kpsi = zeros( ComplexF64, size(psi) )
    for ist = 1:Nstates
        for ip = 1:Npoints
            Kpsi[ip,ist] = psi[ip,ist] / ( 1.0 + G2[ip] )
        end
    end
    return Kpsi
end

# no preconditioning
function Kprec(psi)
  return psi
end
