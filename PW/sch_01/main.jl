include("../common/PWGrid_v01.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("gen_dr.jl")
include("init_pot_harm_3d.jl")
include("apply_K.jl")
include("apply_Vpot.jl")
include("apply_H.jl")
include("calc_rho.jl")
include("gradE.jl")
include("calc_Etot.jl")
include("schsolve_Emin_sd.jl")
include("schsolve_Emin_cg.jl")
include("Kprec.jl")

function test_main( ns1::Int,ns2::Int,ns3::Int )
  #
  Ns = [ns1,ns2,ns3]
  const LatVecs = diagm( [6.0, 6.0, 6.0] )
  #
  pw = PWGrid( Ns, LatVecs )
  #
  const Npoints = pw.Npoints
  const Ω  = pw.Ω
  const r  = pw.r
  const G  = pw.G
  const G2 = pw.G2
  #
  # Generate array of distances
  #
  center = sum(LatVecs,2)/2
  dr = gen_dr( r, center )
  #
  # Setup potential
  #
  Vpot = init_pot_harm_3d( pw, dr )
  print("sum(Vpot)*Ω/Npoints = $(sum(Vpot)*Ω/Npoints)\n");
  #
  const Nstates = 4
  srand(2222)
  psi  = randn(Npoints,Nstates) + im*randn(Npoints,Nstates)
  psi = ortho_gram_schmidt(psi)
  #
  psi, Etot = schsolve_Emin_sd( pw, Vpot, psi, NiterMax=10 )
  psi, Etot = schsolve_Emin_cg( pw, Vpot, psi, NiterMax=1000 )
  #
  Y = ortho_gram_schmidt(psi)
  mu = Y' * apply_H( pw, Vpot, Y )
  evals, evecs = eig(mu)
  Psi = Y*evecs
  for st = 1:Nstates
    @printf("=== State # %d, Energy = %f ===\n", st, real(evals[st]))
  end
end

@time test_main( 30,30,30 )
