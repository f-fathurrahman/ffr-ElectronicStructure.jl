include("../common/PWGrid_v01.jl")
include("../common/gen_lattice.jl")
include("gen_dr.jl")

function test_main(Ns; verbose=false)

  #LatVecs = 6.0*diagm(ones(3))
  #LatVecs = gen_lattice_hexagonal( 7.0, ca=1.5 )
  LatVecs = gen_lattice_fcc( 4.0 )
  #LatVecs = gen_lattice_bcc( 4.0 )
  println(LatVecs)

  pw = PWGrid( Ns, LatVecs )

  dr = gen_dr( pw.R, LatVecs )

  lenx = 0.5*sum(LatVecs[:,1])
  leny = 0.5*sum(LatVecs[:,2])
  lenz = 0.5*sum(LatVecs[:,3])
  MaxDist = sqrt( lenx^2 + leny^2 + lenz^2 )
  @printf("VecLen: (%18.10f,%18.10f,%18.10f)\n", lenx, leny, lenz)
  @printf("Max dist, dr: %18.10f , %18.10f\n", MaxDist, maximum(dr))
end

test_main([4,4,4], verbose=true)
#test_main([75,75,75])
