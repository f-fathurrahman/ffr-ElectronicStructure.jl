include("PWGrid_v02.jl")
include("write_xsf.jl")

function latvec_hexagonal(a; ca=1.0)
  LL = zeros(3,3)
  LL[1,:] = [1.0, 0.0, 0.0]
  LL[2,:] = [cos(pi/3.0), sin(pi/3.0), 0.0]
  LL[3,:] = [0.0, 0.0, ca]
  return a*LL
end

function test_main()
  const ecutwfc = 2 * 0.5  # in Ha
  #LatVecs = 16.0*diagm(ones(3))
  LatVecs = latvec_hexagonal( 16.0, ca=2.0 )

  pw = PWGrid( ecutwfc, LatVecs )

  Ns = pw.Ns

  atpos = pw.R
  write_xsf("R_grid.xsf", LatVecs, atpos)

  scal = 2.0
  Rec = pw.RecVecs*scal
  atpos = pw.gvectors.G*scal
  write_xsf("G_grid.xsf", Rec, atpos, molecule=true)
end

test_main()
