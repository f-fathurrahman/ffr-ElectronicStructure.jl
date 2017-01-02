include("../common/PWGrid_v01.jl")

function write_XSF( filnam, LL, atpos; molecule=false )
  #
  f = open(filnam, "w")
  Natoms = size(atpos)[2]
  #
  if molecule
    @printf(f, "MOLECULE\n")
  else
    @printf(f, "CRYSTAL\n")
  end
  @printf(f, "PRIMVEC\n")
  @printf(f, "%18.10f %18.10f %18.10f\n", LL[1,1], LL[1,2], LL[1,3])
  @printf(f, "%18.10f %18.10f %18.10f\n", LL[2,1], LL[2,2], LL[2,3])
  @printf(f, "%18.10f %18.10f %18.10f\n", LL[3,1], LL[3,2], LL[3,3])
  @printf(f, "PRIMCOORD\n")
  @printf(f, "%8d %8d\n", Natoms, 1)
  for ia = 1:Natoms
    @printf(f, "X  %18.10f %18.10f %18.10f\n", atpos[1,ia], atpos[2,ia], atpos[3,ia])
  end
  close(f)
end


function latvec_hexagonal(a; ca=1.0)
  LL = zeros(3,3)
  LL[1,:] = [1.0, 0.0, 0.0]
  LL[2,:] = [cos(pi/3.0), sin(pi/3.0), 0.0]
  LL[3,:] = [0.0, 0.0, ca]
  return a*LL
end

function latvec_fcc(a)
  LL = zeros(3,3)
  LL[1,:] = [-1.0, 0.0, 1.0]
  LL[2,:] = [ 0.0, 1.0, 1.0]
  LL[3,:] = [-1.0, 1.0, 0.0]
  return 0.5*a*LL
end


function latvec_bcc(a)
  LL = zeros(3,3)
  LL[1,:] = [ 1.0,  1.0, 1.0]
  LL[2,:] = [-1.0,  1.0, 1.0]
  LL[3,:] = [-1.0, -1.0, 1.0]
  return 0.5*a*LL
end


function test_main()
  Ns = [15, 15, 15]

  #LL = 16.0*diagm([1.0, 1.0, 1.0])
  LL = latvec_hexagonal(16.0, ca=2.0)
  #LL = latvec_fcc(16.0)
  #LL = latvec_bcc(16.0)

  pwgrid = PWGrid( Ns, LL )
  atpos = pwgrid.r

  write_XSF("R_grid.xsf", LL, atpos)

  Rec = pwgrid.RecVecs*Ns[1]/2.0
  atpos = pwgrid.G
  write_XSF("G_grid.xsf", LL, atpos, molecule=true)
end


test_main()
