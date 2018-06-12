using Printf
using LinearAlgebra
using Random

include("../common/PWGrid_v02.jl")

include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")
include("../common/gen_lattice_pwscf.jl")

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

function test_main( Ns; solution_method="diag_lobpcg" )

    LatVecs = gen_lattice_sc(6.0)
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

    # Generate array of distances
    center = sum(LatVecs,dims=2)/2
    dr = gen_dr( r, center )

    # Setup potential
    Vpot = init_pot_harm_3d( pw, dr )
    println("sum(Vpot)*Ω/Npoints = ", sum(Vpot)*Ω/Npoints)

    Nstates = 4
    srand(2222)
    psi  = rand(ComplexF64,Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)

    if solution_method == "Emin"
        #psi, Etot = schsolve_Emin_sd( pw, Vpot, psi, NiterMax=10 )
        psi, Etot = schsolve_Emin_cg( pw, Vpot, psi, NiterMax=1000 )

        Y = ortho_gram_schmidt(psi)
        mu = Y' * op_H( pw, Vpot, Y )
        evals, evecs = eigen(mu)
        Psi = Y*evecs

    elseif solution_method == "lobpcg"
        evals, psi = diag_lobpcg( pw, Vpot, psi, verbose=true, tol_avg=1e-10 )

    elseif solution_method == "davidson"
        evals, psi = diag_davidson( pw, Vpot, psi, verbose=true, tol_avg=1e-10, NiterMax=20 )
        evals, psi = diag_davidson( pw, Vpot, psi, verbose=true, tol_avg=1e-10, NiterMax=100 )

    else
        println("Unknown solution_method = ", solution_method)
        exit()

    end

    for ist = 1:Nstates
        @printf("State # %d, Energy = %18.10f\n", ist, real(evals[ist]))
    end

end

@time test_main( [30, 30, 30], solution_method="Emin" )
#@time test_main( [30, 30, 30], solution_method="lobpcg" )
#@time test_main( [30, 30, 30], solution_method="davidson" )