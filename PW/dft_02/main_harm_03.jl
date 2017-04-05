include("../common/PWGrid_v03.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("EnergiesT.jl")
include("PotentialsT.jl")
include("gen_dr.jl")
include("init_pot_harm_3d.jl")
include("op_K.jl")
include("op_V_loc.jl")
include("op_H.jl")
include("calc_rho.jl")
include("calc_grad.jl")
include("calc_Energies.jl")
include("KS_solve_Emin_sd.jl")
include("KS_solve_Emin_cg.jl")
include("Poisson_solve.jl")
include("LDA_VWN.jl")
include("Kprec.jl")

function test_main( ecutwfc_Ry::Float64 )

  const LatVecs = 6.0*diagm( ones(3) )

  pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs )

  const Ω  = pw.Ω
  const r  = pw.r
  const G  = pw.gvec.G
  const G2 = pw.gvec.G2
  const Ns = pw.Ns
  const Npoints = prod(Ns)
  const Ngwx = pw.gvecw.Ngwx

  @printf("ecutwfc (Ha) = %f\n", ecutwfc_Ry*0.5)
  @printf("Ns   = (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])
  @printf("Ngwx = %d\n", Ngwx)

  const actual = Npoints/Ngwx
  const theor = 1/(4*pi*0.25^3/3)
  @printf("Compression: actual, theor: %f , %f\n", actual, theor)

  #
  # Generate array of distances
  #
  center = sum(LatVecs,2)/2
  dr = gen_dr( r, center )
  #
  # Setup potential
  #
  V_ionic = init_pot_harm_3d( pw, dr )
  print("sum(Vpot)*Ω/Npoints = $(sum(V_ionic)*Ω/Npoints)\n");
  #
  const Nstates = 4
  Focc = 2.0*ones(Nstates)
  #
  #psi, Energies, Potentials = KS_solve_Emin_sd( pw, V_ionic, Focc, Nstates, NiterMax=10 )
  #psi, Energies, Potentials = KS_solve_Emin_cg( pw, V_ionic, Focc, Nstates,
  #                            NiterMax=1000, Potentials0=Potentials, psi0=psi )

  psi, Energies, Potentials = KS_solve_Emin_cg( pw, V_ionic, Focc, Nstates, NiterMax=1000 )

  #
  print_Energies(Energies)

  Y = ortho_gram_schmidt(psi)
  mu = Y' * op_H( pw, Potentials, Y )
  evals, evecs = eig(mu)
  Psi = Y*evecs
  for st = 1:Nstates
    @printf("=== State # %d, Energy = %f ===\n", st, real(evals[st]))
  end
end

@time test_main( 40.0 )
