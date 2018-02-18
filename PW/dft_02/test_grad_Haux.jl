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
include("calc_grad_v2.jl")
include("calc_rho.jl")
include("Poisson_solve.jl")
include("LDA_VWN.jl")

include("calc_grad_Haux.jl")
include("calc_Focc.jl")
include("smear_FD.jl")

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
    const Nstates = 5
    const Nelectrons = 4.0
    Focc = zeros(Nstates)
    Focc[1:2] = [2.0, 2.0]

    srand(1234)
    psi = randn(Ngwx,Nstates) + im*randn(Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)

    Potentials = PotentialsT( V_ionic, zeros(Npoints), zeros(Npoints) )
    rho = calc_rho( pw, Focc, psi )
    Potentials.Hartree = real( G_to_R( Ns, Poisson_solve(pw, rho) ) )
    Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )

    # Get eigenvalues
    mu = psi' * op_H(pw, Potentials, psi)
    evals = sort(real(eigvals(mu)))
    println(evals)

    kT = 0.01
    Focc, E_fermi = calc_Focc(evals, Nelectrons, kT, is_spinpol=false)
    
    println("Focc = ", Focc)

    g_Haux = calc_grad_Haux( pw, Potentials, Focc, evals, psi, 0.01 )

    println("sum g_Haux = ", sum(g_Haux))

end

test_main( [30,30,30] )
