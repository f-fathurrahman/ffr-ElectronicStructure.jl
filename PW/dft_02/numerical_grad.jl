include("../common/PWGrid_v02.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("EnergiesT.jl")
include("PotentialsT.jl")
include("gen_dr.jl")
include("init_pot_harm_3d.jl")
include("op_K.jl")
include("op_V_loc.jl")
include("op_H.jl")
include("calc_grad.jl")
include("calc_rho.jl")
include("Poisson_solve.jl")
include("LDA_VWN.jl")
include("calc_Energies.jl")

function test_main( Ns )

  const LatVecs = 16.0*diagm( ones(3) )

  pw = PWGrid( Ns, LatVecs )

  const Ω  = pw.Ω
  const r  = pw.r
  const G  = pw.gvec.G
  const G2 = pw.gvec.G2
  const Npoints = prod(Ns)
  const Ngwx = pw.gvecw.Ngwx

  #
  # Generate array of distances
  #
  center = sum(LatVecs,2)/2
  dr = gen_dr( r, center )
  #
  # Setup potential
  #
  V_ionic = init_pot_harm_3d( pw, dr )
  #
  const Nstates = 4
  Focc = 2.0*ones(Nstates)

  # Initialize random wave function and orthonormalize in G-space
  srand(1234)
  psi = randn(Ngwx,Nstates) + im*randn(Ngwx,Nstates)
  psi = ortho_gram_schmidt(psi)

  Potentials = PotentialsT( V_ionic, zeros(Npoints), zeros(Npoints) )
  rho = calc_rho( pw, Focc, psi )
  Potentials.Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
  Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )

  Energies = calc_Energies( pw, Potentials, Focc, psi, 0.0 )
  E0 = Energies.Total
  g0 = calc_grad( pw, Potentials, Focc, psi )

  println("E0 = ", E0)
  println("sum(g0) = ", sum(g0))

  # random direction
  dW = randn(Ngwx,Nstates) + im*randn(Ngwx,Nstates)

  for ii = 1:9
    delta = 10.0^(-ii)

    psiw = psi + delta*dW
    psiw = ortho_gram_schmidt(psiw)
    Potentials = PotentialsT( V_ionic, zeros(Npoints), zeros(Npoints) )
    rho = calc_rho( pw, Focc, psiw )
    Potentials.Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
    Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )

    dE = 2*real(trace(g0'*delta*dW))

    Energies = calc_Energies( pw, Potentials, Focc, psiw, 0.0 )
    E = Energies.Total
    g = calc_grad( pw, Potentials, Focc, psiw )

    println("\ndelta = ", delta)
    # actual change / expected change
    println("Actual / expected = ", (E-E0)/dE)
    # e
    println("Error estimate = ", sqrt(Ngwx)*eps()/abs(dE) )

  end

end

test_main( [30,30,30] )
