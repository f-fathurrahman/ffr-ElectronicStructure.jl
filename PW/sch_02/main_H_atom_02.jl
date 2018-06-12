using Printf
using LinearAlgebra
using Random

include("../common/PWGrid_v02.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")
include("../common/gen_lattice_pwscf.jl")

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

function test_main( Ns )

    LatVecs = gen_lattice_sc(16.0)

    pw = PWGrid( Ns, LatVecs )

    Ω  = pw.Ω
    r  = pw.r
    G  = pw.gvec.G
    G2 = pw.gvec.G2
    Npoints = prod(Ns)
    Ngwx = pw.gvecw.Ngwx

    @printf("Ns   = (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])
    @printf("Ngwx = %d\n", Ngwx)

    actual = Npoints/Ngwx
    theor = 1/(4*pi*0.25^3/3)
    @printf("Compression: actual, theor: %f , %f\n", actual, theor)

    Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )

    Sf = structure_factor( Xpos, G )

    E_nn = calc_ewald( pw, Xpos, Sf )

    Vg = zeros(ComplexF64,Npoints)
    prefactor = -4*pi/Ω
    for ig=2:Npoints
        Vg[ig] = prefactor/G2[ig]
    end
    Vpot = real( G_to_R(Ns, Vg .* Sf) ) * Npoints

    Nstates = 1
    srand(2222)
    psi = rand(ComplexF64,Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)
    #
    #evals, evecs = diag_lobpcg( pw, Vpot, psi, verbose=true, tol_avg=1e-10 )

    psi, Etot = schsolve_Emin_sd( pw, Vpot, psi, NiterMax=10 )
    psi, Etot = schsolve_Emin_cg( pw, Vpot, psi, NiterMax=1000 )
    #
    Y = ortho_gram_schmidt(psi)
    mu = Y' * op_H( pw, Vpot, Y )
    evals, evecs = eigen(mu)
    psi = Y*evecs

    for ist = 1:Nstates
        @printf("State %d, Energy = %18.10f\n", ist, real(evals[ist]))
    end

    @printf("E_nn    = %18.10f\n", E_nn)
    @printf("E total = %18.10f\n", E_nn + Etot)

end

@time test_main( [64, 64, 64] )
@time test_main( [64, 64, 64] )
