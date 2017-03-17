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
include("kssolve_Emin_sd.jl")
include("kssolve_Emin_cg.jl")
include("solve_poisson.jl")
include("LDA_VWN.jl")
include("Kprec.jl")

include("structure_factor.jl")
include("gen_rho.jl")
include("gen_dr.jl")
include("calc_ewald.jl")

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

  Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )

  Sf = structure_factor( Xpos, G )

  E_nn = calc_ewald( pw, Xpos, Sf )

  Vg = zeros(Complex128,Npoints)
  prefactor = -4*pi/Ω
  for ig=2:Npoints
    Vg[ig] = prefactor/G2[ig]
  end
  V_ionic = real( G_to_R(Ns, Vg .* Sf) ) * Npoints

  println("sum(Sf) = $(sum(Sf))")
  println("sum(Vg) = $(sum(Vg))")
  println("sum(V_ionic) = $(sum(V_ionic))")
  for ip = 1:10
    @printf("%8d %18.10f\n", ip, V_ionic[ip])
  end
  @printf("Ω = %f\n", Ω)
  @printf("maximum(Vpot) = %18.10f\n", maximum(V_ionic))
  @printf("minimum(Vpot) = %18.10f\n", minimum(V_ionic))

  #
  const Nstates = 1
  Focc = [1.0]
  psi, Energies, Potentials = kssolve_Emin_cg( pw, V_ionic, Focc, Nstates, NiterMax=1000 )

  #
  Y = ortho_gram_schmidt(psi)
  mu = Y' * op_H( pw, Potentials, Y )
  evals, evecs = eig(mu)
  psi = Y*evecs

  for st = 1:Nstates
    @printf("=== State # %d, Energy = %f ===\n", st, real(evals[st]))
  end

  @printf("E_nn    = %18.10f\n", E_nn)
  @printf("E total = %18.10f\n", E_nn + Energies.Total)

end

@code_native test_main(1.0)
@time test_main( 35.0 )
