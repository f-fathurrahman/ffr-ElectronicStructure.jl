include("../common/PWGrid_v02.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("EnergiesT.jl")
include("PotentialsT.jl")
include("gen_dr.jl")
include("apply_K.jl")
include("apply_V_loc.jl")
include("apply_H.jl")
include("calc_rho.jl")
include("gradE.jl")
include("calc_Energies.jl")
include("kssolve_Emin_sd.jl")
include("kssolve_Emin_cg.jl")
include("solve_poisson.jl")
include("LDA_VWN.jl")
include("Kprec.jl")

include("../common/gen_dr_center.jl")
include("../common/calc_strfact_v1.jl")
include("../common/calc_ewald_v1.jl")

include("diag_lobpcg.jl")
include("kssolve_scf.jl")

function test_main( Ns )

  const a = 5.66/0.52917721 # Lattice constant (converted from angstroms to bohrs)
  const LatVecs = a*diagm(ones(3))

  pw = PWGrid( Ns, LatVecs )

  const Ω  = pw.Ω
  const r  = pw.r
  const G  = pw.gvec.G
  const G2 = pw.gvec.G2
  const Npoints = prod(Ns)
  const Ngwx = pw.gvecw.Ngwx

  @printf("Ns   = (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])
  @printf("Ngwx = %d\n", Ngwx)

  const actual = Npoints/Ngwx
  const theor = 1/(4*pi*0.25^3/3)
  @printf("Compression: actual, theor: %f , %f\n", actual, theor)

  # diamond lattice in cubic cell
  Xpos = zeros(3,8)
  Xpos[:,1] = [0.00 0.00 0.00]
  Xpos[:,2] = [0.25 0.25 0.25]
  Xpos[:,3] = [0.00 0.50 0.50]
  Xpos[:,4] = [0.25 0.75 0.75]
  Xpos[:,5] = [0.50 0.00 0.50]
  Xpos[:,6] = [0.75 0.25 0.75]
  Xpos[:,7] = [0.50 0.50 0.00]
  Xpos[:,8] = [0.75 0.75 0.25]
  Xpos = a*Xpos[:,:]

  for ia = 1:size(Xpos)[2]
    @printf("Ge %18.10f %18.10f %18.10f\n", Xpos[1,ia], Xpos[2,ia], Xpos[3,ia] )
  end

  println(LatVecs)

  Nspecies = 1
  atmsymb = ["Ge"] # unique list of atomic symbols
  atm2species = [1, 1, 1, 1, 1, 1, 1, 1]  # mapping from atom to species
  Zv = [4.0]  # only valence ?

  Sf = calc_strfact( Xpos, Nspecies, atm2species, pw.gvec.G )

  E_nn = calc_ewald( pw, Sf, Xpos, Nspecies, atm2species, Zv )
  @printf("E_nn = %18.10f\n", E_nn)

  Vps = zeros(Complex128, Npoints)

  # Ge pseudopotential
  Z = Zv[1]
  λ = 18.5
  rc = 1.052
  Gm = sqrt(G2)

  Vps[:] = -2*pi*exp(-pi*Gm/λ).*cos(Gm*rc).*(Gm/λ)./(1-exp(-2*pi*Gm/λ))
  for n = 0:4
    Vps[:] = Vps[:] + (-1)^n*exp(-λ*rc*n)./(1+(n*λ./Gm).^2)
  end

  Vps[:] = Vps.*4*pi*Z./Gm.^2*(1+exp(-λ*rc)) - 4*pi*Z./Gm.^2

  n = collect(1:4)
  Vps[1] = 4*pi*Z*(1+exp(-λ*rc))*(rc^2/2+1/λ^2 *
           (pi^2/6 + sum((-1).^n.*exp(-λ*rc*n)./n.^2)))

  Vps = Vps[:]/Ω

  V_ionic = real( G_to_R(Ns, Vps .* Sf) ) * Npoints

  # Need to sum up over Nspecies for more than one species type
  # We simply need reshape because we only have one species type here.
  V_ionic = reshape( V_ionic, (Npoints) )

  const Nstates = 16
  Focc = 2.0*ones(Nstates)

  psi, Energies, Potentials = kssolve_Emin_cg( pw, V_ionic, Focc, Nstates, NiterMax=1000 )
  #psi, Energies, Potentials = kssolve_Emin_cg( pw, V_ionic, Focc, Nstates,
  #                            NiterMax=1000, Potentials0=Potentials, psi0=psi )

  Y = ortho_gram_schmidt(psi)
  mu = Y' * apply_H( pw, Potentials, Y )
  evals, evecs = eig(mu)
  psi = Y*evecs

  #Energies, Potentials, psi, evals = kssolve_scf( pw, V_ionic, Focc, Nstates )

  for st = 1:Nstates
    @printf("=== State # %4d, Energy = %18.10f ===\n", st, real(evals[st]))
  end

  @printf("E_nn    = %18.10f\n", E_nn)
  @printf("E total = %18.10f\n", E_nn + Energies.Total)

  @printf("\nE total (eV per atom)= %18.10f\n", (E_nn + Energies.Total)/8*27.21 )

end

#@time test_main( [48,48,48] )
test_main( [46, 46, 46] )

