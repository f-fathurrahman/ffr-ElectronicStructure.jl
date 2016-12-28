include("../common/PWGrid_v01.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("apply_K.jl")
include("apply_Vpot.jl")
include("apply_H.jl")
include("calc_rho.jl")
include("calc_Etot.jl")
include("Kprec.jl")
include("diag_lobpcg.jl")
include("gradE.jl")
include("schsolve_Emin_sd.jl")
include("schsolve_Emin_cg.jl")
include("structure_factor.jl")

function test_main( ns1::Int,ns2::Int,ns3::Int )

  Ns = [ns1,ns2,ns3]
  const LatVecs = 70.0*diagm( ones(3) )

  pw = PWGrid( Ns, LatVecs )

  const Npoints = pw.Npoints
  const Ω = pw.Ω
  const R = pw.R
  const G = pw.G
  const G2 = pw.G2

  Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )

  Sf = structure_factor( Xpos, G )
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
  psi = randn(Npoints,Nstates) + im*randn(Npoints,Nstates)
  psi = ortho_gram_schmidt(psi)
  #
  #evals, evecs = diag_lobpcg( pw, Vpot, psi, verbose=true, tol_avg=1e-10 )

  psi, Etot = schsolve_Emin_sd( pw, Vpot, psi, NiterMax=10 )
  psi, Etot = schsolve_Emin_cg( pw, Vpot, psi, NiterMax=1000 )
  #
  Y = ortho_gram_schmidt(psi)
  mu = Y' * apply_H( pw, Vpot, Y )
  evals, evecs = eig(mu)
  psi = Y*evecs

  for st = 1:Nstates
    @printf("=== State # %d, Energy = %f ===\n", st, real(evals[st]))
  end

end

@time test_main( 100,100,100 )
