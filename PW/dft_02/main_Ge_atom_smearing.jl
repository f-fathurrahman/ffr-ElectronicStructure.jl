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
include("KS_solve_Emin_sd.jl")
include("KS_solve_Emin_cg.jl")
include("Poisson_solve.jl")
include("LDA_VWN.jl")
include("Kprec.jl")

include("../common/gen_dr_center.jl")
include("../common/calc_strfact_v2.jl")
include("../common/calc_ewald_v2.jl")

include("diag_lobpcg.jl")
include("KS_solve_SCF_smearing.jl")

include("getocc.jl")
include("fermidirac.jl")
include("getEntropy.jl")

include("andersonmix.jl")

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

    @printf("Ns = (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])
    @printf("Ngwx = %d\n", Ngwx)

    const actual = Npoints/Ngwx
    const theor = 1/(4*pi*0.25^3/3)
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

    Vps = zeros(Complex128, Npoints)

    # Ge pseudopotential
    Z = Zv[1]
    λ = 18.5
    rc = 1.052
    Gm = sqrt.(G2)

    Vps[:] = -2*pi*exp.(-pi*Gm/λ).*cos.(Gm*rc).*(Gm/λ)./(1-exp.(-2*pi*Gm/λ))
    for n = 0:4
        Vps[:] = Vps[:] + (-1)^n*exp(-λ*rc*n)./(1+(n*λ./Gm).^2)
    end

    Vps[:] = Vps.*4*pi*Z./Gm.^2*(1+exp(-λ*rc)) - 4*pi*Z./Gm.^2

    n = collect(1:4)
    Vps[1] = 4*pi*Z*(1+exp(-λ*rc))*(rc^2/2+1/λ^2 *
                     (pi^2/6 + sum((-1).^n.*exp.(-λ*rc*n)./n.^2)))

    Vps = Vps[:]/Ω

    V_ionic = real( G_to_R(Ns, Vps .* Sf) ) * Npoints

    # needed to sum up over Nspecies for more than one species
    V_ionic = reshape( V_ionic, (Npoints) )

    const Nstates = 6
    # Initial Focc
    Focc = zeros(Nstates)
    Focc[1] = 2.0
    Focc[2] = 2.0
    Nocc = 4

    Energies, Potentials, psi, evals = KS_solve_SCF_smearing( pw, V_ionic, Focc, Nstates, Nocc, β=0.5 )

    for st = 1:Nstates
        @printf("=== State # %4d, Energy = %18.10f ===\n", st, real(evals[st]))
    end

    @printf("E_nn = %18.10f\n", E_nn)
    @printf("E total = %18.10f\n", E_nn + Energies.Total)

end

@time test_main( [24,24,24] )
