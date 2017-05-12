include("../common/PWGrid_v03.jl")

include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("gen_dr.jl")
include("init_pot_harm_3d.jl")
include("op_K.jl")
include("op_Vpot.jl")
include("op_H.jl")
include("calc_rho.jl")
include("calc_grad.jl")
include("calc_Etot.jl")
include("schsolve_Emin_sd.jl")
include("schsolve_Emin_cg.jl")
include("Kprec.jl")
include("diag_lobpcg.jl")
include("diag_davidson.jl")

function test_main( ecutwfc_Ry; solution_method="diag_lobpcg" )

  const LatVecs = 6.0*diagm( ones(3) )

  pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs )

  const Ω  = pw.Ω
  const r  = pw.r
  const G  = pw.gvec.G
  const G2 = pw.gvec.G2
  const Ns = pw.Ns
  const Npoints = prod(Ns)
  const Ngwx = pw.gvecw.Ngwx

  @printf("Ns   = (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])
  @printf("Ngwx = %d\n", Ngwx)

  const actual = Npoints/Ngwx
  const theor = 1/(4*pi*0.25^3/3)
  @printf("Compression: actual, theor: %f , %f\n", actual, theor)

  # Generate array of distances
  center = 6.0*ones(3)/2
  dr = gen_dr( r, center )

  # Setup potential
  Vpot = init_pot_harm_3d( pw, dr )
  print("sum(Vpot)*Ω/Npoints = $(sum(Vpot)*Ω/Npoints)\n")

  #
  const Nstates = 4
  srand(2222)
  psi  = randn(Ngwx,Nstates) + im*randn(Ngwx,Nstates)
  psi = ortho_gram_schmidt(psi)
  #psi = zeros(Complex128,Ngwx,Nstates)
  #for ist = 1:Nstates
  #  psi[ist,ist] = 1.0
  #end

  #Kpsi = op_K( pw, psi )
  #println("sum(Kpsi) = ", sum(Kpsi))

  #Vpsi = op_Vpot( pw, Vpot, psi )
  #println("sum(Vpsi) = ", sum(Vpsi))

  #Hpsi = op_H( pw, Vpot, psi )
  #println("sum(Hpsi) = ", sum(Hpsi))

  #grad = calc_grad( pw, Vpot, psi )
  #println("sum(grad) = ", sum(grad))

  #Etot = calc_Etot( pw, Vpot, psi )
  #println("Etot = ", Etot)

  psi, Etot = schsolve_Emin_sd( pw, Vpot, psi, NiterMax=10 )
  psi, Etot = schsolve_Emin_cg( pw, Vpot, psi, NiterMax=100 )

end


test_main( 50.0 )
