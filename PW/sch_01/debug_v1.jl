#FFTW.set_num_threads(2);
#blas_set_num_threads(2);

include("PWGrid.jl")
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
include("gradE.jl")
include("get_Etot.jl")
include("Etot_minimize_sd.jl")
include("Etot_minimize_cg.jl")

function test_main( ns1::Int,ns2::Int,ns3::Int )
  #
  Ns = [ns1,ns2,ns3]
  const LatVecs = diagm( [6.0, 6.0, 6.0] )
  #
  pw_grid = PWGrid( Ns, LatVecs )
  #
  const Npoints = pw_grid.Npoints
  const Ω = pw_grid.Ω
  const R = pw_grid.R
  const G = pw_grid.G
  const G2 = pw_grid.G2
  #
  # Generate array of distances
  #
  center = sum(LatVecs,2)/2
  dr = gen_dr( R, center )
  #
  # Setup potential
  #
  Vpot = gen_Vpot( pw_grid, dr )
  print("sum(Vpot)*Ω/Npoints = $(sum(Vpot)*Ω/Npoints)\n");
  #
  const Nstates = 4
  srand(2222)
  psi  = randn(Npoints,Nstates) + im*randn(Npoints,Nstates)
  #for is = 1:Nstates
  #  psi[is,is] = 1.0
  #end
  # orthonormalization, according to: dot(psi,psi)*Omega = 1
  #psi = scale( 1.0/sqrt(Ω), psi*inv(sqrtm(psi'*psi)) )
  #ortho_gram_schmidt!( Nstates, psi ); scale!( 1./sqrt(Ω), psi )
  psi = ortho_gram_schmidt(psi)/sqrt(Ω)
  #test_ortho( 1.0, Nstates, psi )
  #Hpsi = apply_H( pw_grid, psi, Vpot )
  #
  #print("Testing energy:\n")
  #Etot = get_Etot( pw_grid, psi, Vpot )
  #Etot_orig = get_Etot_orig( pw_grid, psi, Vpot )
  #print("Testing gradient:\n")
  #Kpsi = apply_K( pw_grid, psi )
  #Vpsi = apply_Vpot( pw_grid, psi, Vpot )
  #Hpsi = Kpsi + Vpsi
  #print("sum(Kpsi) = $(sum(Kpsi))\n")
  #print("sum(Vpsi) = $(sum(Vpsi))\n")
  #print("sum(Hpsi) = $(sum(Hpsi))\n")
  #∇E = gradE( pw_grid, psi, Vpot )
  #print("sum(∇E) = $(sum(∇E))\n")
  #∇E_v2 = gradE_v2( pw_grid, psi, Vpot )
  #print("sum(∇E_v2) = $(sum(∇E_v2))\n")
  ##print("ratio = $(∇E_v2[2]/∇E[2])\n")
  #∇E_norm = gradE_orthonorm( pw_grid, psi, Vpot )
  #print("sum(∇E_norm) = $(sum(∇E_norm))\n")
  #exit()
  #
  psi, Etot = Etot_minimize_sd( pw_grid, psi, Vpot, 10 )
  psi, Etot = Etot_minimize_cg( pw_grid, psi, Vpot, 1000 )
  #
  Y = copy(psi)
  ortho_gram_schmidt!( Nstates, Y ); scale!( 1./sqrt(Ω), Y )
  mu = Y' * apply_H( pw_grid, Y, Vpot )
  evals, evecs = eig(mu)
  Psi = Y*evecs
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
