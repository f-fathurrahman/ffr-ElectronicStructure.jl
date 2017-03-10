include("../common/PWGrid_v03.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("op_K.jl")
include("op_Vpot.jl")
include("op_H.jl")
include("calc_rho.jl")
include("calc_Etot.jl")
include("Kprec.jl")
include("diag_lobpcg.jl")
include("calc_grad.jl")
include("schsolve_Emin_sd.jl")
include("schsolve_Emin_cg.jl")
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
  Vpot = real( G_to_R(Ns, Vg .* Sf) ) * Npoints

  println("sum(Sf) = $(sum(Sf))")
  println("sum(Vg) = $(sum(Vg))")
  println("sum(Vpot) = $(sum(Vpot))")
  for ip = 1:10
    @printf("%8d %18.10f\n", ip, Vpot[ip])
  end
  @printf("Ω = %f\n", Ω)
  @printf("maximum(Vpot) = %18.10f\n", maximum(Vpot))
  @printf("minimum(Vpot) = %18.10f\n", minimum(Vpot))

  #
  const Nstates = 1
  srand(2222)
  psi = randn(Ngwx,Nstates) + im*randn(Ngwx,Nstates)
  psi = ortho_gram_schmidt(psi)
  #
  #evals, evecs = diag_lobpcg( pw, Vpot, psi, verbose=true, tol_avg=1e-10 )

  psi, Etot = schsolve_Emin_sd( pw, Vpot, psi, NiterMax=10 )
  psi, Etot = schsolve_Emin_cg( pw, Vpot, psi, NiterMax=1000 )
  #
  Y = ortho_gram_schmidt(psi)
  mu = Y' * op_H( pw, Vpot, Y )
  evals, evecs = eig(mu)
  psi = Y*evecs

  for st = 1:Nstates
    @printf("=== State # %d, Energy = %f ===\n", st, real(evals[st]))
  end

  @printf("E_nn    = %18.10f\n", E_nn)
  @printf("E total = %18.10f\n", E_nn + Etot)

end

@time test_main( 60.0 )
