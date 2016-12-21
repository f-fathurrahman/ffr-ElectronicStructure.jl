#FFTW.set_num_threads(2);
#blas_set_num_threads(2);

include("EnergiesT.jl")
include("PotentialsT.jl")

include("PWGrid.jl")
include("R_to_G.jl")
include("G_to_R.jl")
include("calc_dr.jl")
include("init_V_harm_3d.jl")
include("apply_K.jl")
include("apply_V_loc.jl")
include("apply_H.jl")
include("ortho_gram_schmidt.jl")
include("test_ortho.jl")
include("calc_rho.jl")
include("solve_poisson.jl")
include("LDA_VWN.jl")
include("calc_grad.jl")
include("calc_Energies.jl")
include("E_min_sd.jl")
include("E_min_pcg.jl")

function test_main( ns1::Int,ns2::Int,ns3::Int )
  #
  Ns = [ns1,ns2,ns3]
  const LatVecs = diagm( [6.0, 6.0, 6.0] )
  #
  PW = PWGrid( Ns, LatVecs )
  #
  const Npoints = PW.Npoints
  const Ω = PW.Ω
  const R = PW.R
  const G = PW.G
  const G2 = PW.G2
  # Generate array of distances
  center = sum(LatVecs,2)/2
  dr = calc_dr( R, center )
  # Setup potential
  V_ionic = init_V_harm_3d( PW, dr )

  Potentials = PotentialsT( V_ionic, zeros(Npoints), zeros(Npoints) )

  Energies = EnergiesT( 0.0, 0.0, 0.0, 0.0, 0.0 )

  #print("sum(V_ionic)*Ω/Npoints = $(sum(V_ionic)*Ω/Npoints)\n");
  #
  const Nstates = 4
  Focc = 2.0*ones(Nstates)
  #
  srand(2222)
  psi  = randn(Npoints,Nstates) + im*randn(Npoints,Nstates)
  #
  # orthonormalization, according to: dot(psiG,psiG) = 1
  #
  #psi = psi*inv(sqrtm(psi'*psi))
  #ortho_gram_schmidt!( Nstates, psi )
  psi = ortho_gram_schmidt(psi)
  #

  Etot, psi = E_min_sd!( PW, Energies, Potentials, Focc, psi, verbose=true, NiterMax=10 )
  Etot, psi = E_min_pcg!( PW, Energies, Potentials, Focc, psi, verbose=true, NiterMax=200 )

  #HH = psi' * apply_H( PW, Potentials, Focc, psi )
  # Alternatively: recalculate local potentials
  rho = calc_rho( PW, Focc, psi )
  Potentials.Hartree = real( G_to_R( PW.Ns, solve_poisson( PW, rho ) ) )
  Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )
  #
  HH = psi' * apply_H( PW, Potentials.Hartree + Potentials.Ionic + Potentials.XC, psi )
  evals, evecs = eig(HH)
  evecs = psi*evecs
  for st = 1:Nstates
    @printf("=== State # %d, Energy = %f ===\n", st, real(evals[st]))
  end
end

#@code_native test_main( 2,2,2 )
@time test_main( 30,30,30 )
