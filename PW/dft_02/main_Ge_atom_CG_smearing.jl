using Printf
using LinearAlgebra
using Random

include("../common/PWGrid_v02.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("EnergiesT.jl")
include("PotentialsT.jl")
include("gen_dr.jl")
include("op_K.jl")
include("op_V_loc.jl")
include("op_H.jl")
include("calc_rho.jl")
include("calc_grad.jl")
include("calc_Energies.jl")
include("KS_solve_Emin_CG_smearing.jl")
include("Poisson_solve.jl")
include("LDA_VWN.jl")
include("Kprec.jl")

include("../common/gen_dr_center.jl")
include("../common/calc_strfact_v2.jl")
include("../common/calc_ewald_v2.jl")

include("calc_Focc.jl")
include("smear_FD.jl")
include("calc_entropy.jl")
include("sum_upto_E_fermi.jl")
include("calc_grad_Haux.jl")
include("andersonmix.jl")

include("printMatrix.jl")

function test_main( Ns )

    Random.seed!(1234)

    a = 5.66/0.52917721 # Lattice constant (converted from angstroms to bohrs)
    LatVecs = a*diagm(0 => ones(3))

    pw = PWGrid( Ns, LatVecs )

    Ω  = pw.Ω
    r  = pw.r
    G  = pw.gvec.G
    G2 = pw.gvec.G2
    Npoints = prod(Ns)
    Ngwx = pw.gvecw.Ngwx

    @printf("Ns = (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])
    @printf("Ngwx = %d\n", Ngwx)

    actual = Npoints/Ngwx
    theor = 1/(4*pi*0.25^3/3)
    @printf("Compression: actual, theor: %f , %f\n", actual, theor)

    # diamond lattice in cubic cell
    Xpos = zeros(3,1)
    Xpos[:,1] = [0.0, 0.0, 0.0]

    Nspecies = 1
    atmsymb = ["Ge"] # unique list of atomic symbols
    atm2species = [1]    # mapping from atom to species
    Zv = [4.0]    # only valence ?

    Sf = calc_strfact( Xpos, Nspecies, atm2species, pw.gvec.G )

    E_nn = calc_ewald( pw, Sf, Xpos, Nspecies, atm2species, Zv )
    @printf("E_nn = %18.10f\n", E_nn)

    Vps = zeros(ComplexF64, Npoints)

    # Ge pseudopotential
    Z = Zv[1]
    λ = 18.5
    rc = 1.052

    Vps = zeros(ComplexF64,Npoints)
    for i = 2:Npoints
        Gm = sqrt(G2[i])
        Vps[i] = -2*pi*exp(-pi*Gm/λ)*cos(Gm*rc)*(Gm/λ)/(1 - exp(-2*pi*Gm/λ))
        for n = 0:4
            Vps[i] = Vps[i] + (-1.0)^n*exp(-λ*rc*n)/(1 + (n*λ/Gm)^2)
        end
        Vps[i] = Vps[i]*4*pi*Z/Gm^2*(1+exp(-λ*rc)) - 4*pi*Z/Gm^2
    end

    n = collect(1:4)
    Vps[1] = 4*pi*Z*(1+exp(-λ*rc))*(rc^2/2+1/λ^2 *
                     (pi^2/6 + sum((-1).^n.*exp.(-λ*rc*n)./n.^2)))

    Vps[:] = Vps[:]/Ω

    V_ionic = real( G_to_R(Ns, Vps .* Sf) ) * Npoints

    # needed to sum up over Nspecies for more than one species
    V_ionic = reshape( V_ionic, (Npoints) )

    Nstates = 7
    # Initial Focc
    Focc = zeros(Nstates)
    Focc[1] = 2.0
    Focc[2] = 1.0
    Focc[3] = 1.0
    Nelectrons = 4

    Energies, Potentials, psi, evals = KS_solve_Emin_CG_smearing( pw, V_ionic, Focc, Nstates, E_NN=E_nn )

end

test_main( [24,24,24] )
