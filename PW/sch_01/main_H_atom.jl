include("../common/PWGrid_v01.jl")
include("R_to_G.jl")
include("G_to_R.jl")
include("gen_dr.jl")
include("gen_Vpot.jl")
include("apply_K.jl")
include("apply_Vpot.jl")
include("apply_H.jl")
include("ortho_gram_schmidt.jl")
include("test_ortho.jl")
include("get_rho.jl")
include("get_Etot.jl")
include("Kprec.jl")
include("diag_lobpcg.jl")

function test_main( ns1::Int,ns2::Int,ns3::Int )
  #
  Ns = [ns1,ns2,ns3]
  const LatVecs = 10.0*diagm( ones(3) )
  #
  pw_grid = PWGrid( Ns, LatVecs )
  #
  const Npoints = pw_grid.Npoints
  const Ω = pw_grid.Ω
  const R = pw_grid.R
  const G = pw_grid.G
  const G2 = pw_grid.G2

  Vg = zeros(Complex128,Npoints)
  prefactor = -4*pi/Ω
  for ig=2:Npoints
    Vg[ig] = prefactor/G2[ig]
  end
  Vpot = real( G_to_R(NN,Vg) ) * Npoints

  print("sum(Vpot)*Ω/Npoints = $(sum(Vpot)*Ω/Npoints)\n");
  #
  const Nstates = 1
  srand(2222)
  psi = randn(Npoints,Nstates) + im*randn(Npoints,Nstates)
  psi = ortho_gram_schmidt(psi)
  #
  evals, evecs = diag_lobpcg( pw_grid, Vpot, psi, verbose=true, tol_avg=1e-10 )

  for st = 1:Nstates
    @printf("=== State # %d, Energy = %f ===\n", st, real(evals[st]))
  end
end

@time test_main( 30,30,30 )
