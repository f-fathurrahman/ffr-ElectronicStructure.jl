include("../common/PWGrid_v01.jl")
include("../common/ortho_gram_schmidt.jl")

include("R_to_G.jl")
include("G_to_R.jl")

include("gen_dr.jl")
include("init_pot_harm_3d.jl")
include("apply_K.jl")
include("apply_Vpot.jl")
include("apply_H.jl")

include("calc_rho.jl")
include("calc_Etot.jl")
include("Kprec.jl")
include("diag_lobpcg.jl")

function test_main( ns1::Int,ns2::Int,ns3::Int )
  #
  Ns = [ns1,ns2,ns3]
  const LatVecs = diagm( [6.0, 6.0, 6.0] )
  #
  pw = PWGrid( Ns, LatVecs )
  #
  const Npoints = pw.Npoints
  const Ω = pw.Ω
  const R = pw.R
  const G = pw.G
  const G2 = pw.G2
  #
  # Generate array of distances
  #
  center = sum(LatVecs,2)/2
  dr = gen_dr( R, center )
  #
  # Setup potential
  #
  Vpot = init_pot_harm_3d( pw, dr )
  print("sum(Vpot)*Ω/Npoints = $(sum(Vpot)*Ω/Npoints)\n");
  #
  const Nstates = 4
  srand(2222)
  psi = randn(Npoints,Nstates) + im*randn(Npoints,Nstates)
  #
  evals, evecs = diag_lobpcg( pw, Vpot, psi, verbose=true, tol_avg=1e-10 )

  for st = 1:Nstates
    @printf("=== State # %d, Energy = %f ===\n", st, real(evals[st]))
  end
end

#@code_native test_main( 2,2,2 )
@time test_main( 30,30,30 )


#test_main()
#Profile.clear()
#@profile test_main()
#
#r = Profile.retrieve();
#f = open("Profile.bin", "w")
#serialize(f, r)
#close(f)
