include("../common/PWGrid_v03.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("EnergiesT.jl")
include("PotentialsT.jl")
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

include("../common/gen_dr_center.jl")
include("../common/calc_strfact_v2.jl")
include("../common/calc_ewald_v2.jl")

include("diag_lobpcg.jl")
include("KS_solve_scf.jl")d

function test_main( ecutwfc_Ry::Float64 )

  const LatVecs = 16.0*diagm( ones(3) )

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

  # C atom
  const Zatm = 6.0;
  const Nstates = 3
  Focc = [2.0, 2.0, 2.0]

  Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )

  Sf = calc_strfact( Xpos, 1, [1], pw.gvec.G )
  E_nn = calc_ewald( pw, Sf, Xpos, 1, [1], [Zatm] )

  Vg = zeros(Complex128,Npoints)
  prefactor = -4*pi*Zatm/Ω
  for ig=2:Npoints
    Vg[ig] = prefactor/G2[ig]
  end
  V_ionic = real( G_to_R(Ns, Vg .* Sf) ) * Npoints

  # Need to sum up over Nspecies for more than one species type
  # We simply need reshape because we only have one species type here.
  V_ionic = reshape( V_ionic, (Npoints) )

  #psi, Energies, Potentials = KS_solve_Emin_cg( pw, V_ionic, Focc, Nstates, NiterMax=1000 )

  #
  #Y = ortho_gram_schmidt(psi)
  #mu = Y' * op_H( pw, Potentials, Y )
  #evals, evecs = eig(mu)
  #psi = Y*evecs

  @printf("Solution by self-consistent field method\n")
  Energies, Potentials, psi, evals = KS_solve_scf( pw, V_ionic, Focc, Nstates )

  #for st = 1:Nstates
  #  @printf("=== State # %d, Energy = %f ===\n", st, real(evals[st]))
  #end

  @printf("E_nn    = %18.10f\n", E_nn)
  @printf("E total = %18.10f\n", E_nn + Energies.Total)

end

@time test_main( 35.0 )
