using Printf
using LinearAlgebra
using Random

include("../common/PWGrid_v01.jl")
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
include("Sch_solve_Emin_sd.jl")
include("Sch_solve_Emin_cg.jl")
include("structure_factor.jl")

include("gen_rho.jl")
include("gen_dr.jl")
include("calc_ewald.jl")

function test_main( ns1::Int64,ns2::Int64,ns3::Int64 )

    Ns = [ns1,ns2,ns3]
    LatVecs = gen_lattice_sc(16.0)

    pw = PWGrid( Ns, LatVecs )

    Npoints = pw.Npoints
    Ω = pw.Ω
    r = pw.r
    G = pw.G
    G2 = pw.G2

    Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )

    Sf = structure_factor( Xpos, G )

    E_nn = calc_ewald( pw, Xpos, Sf )

    Vg = zeros(ComplexF64,Npoints)
    prefactor = -4*pi/Ω
    for ig=2:Npoints
        Vg[ig] = prefactor/G2[ig]
    end
    Vpot = real( G_to_R(Ns, Vg .* Sf) ) * Npoints

    println("sum(Sf) = $(sum(Sf))")
    println("sum(Vg) = $(sum(Vg))")
    println("sum(Vpot) = $(sum(Vpot))")
    for ip = 1:10
        @printf("%8d %18.10f\n", ip, Vpot[ip])
    end
    @printf("Ω = %f\n", Ω)
    @printf("maximum(Vpot) = %18.10f\n", maximum(Vpot))
    @printf("minimum(Vpot) = %18.10f\n", minimum(Vpot))

    #
    Nstates = 1
    srand(2222)
    psi = randn(Npoints,Nstates) + im*randn(Npoints,Nstates)
    psi = ortho_gram_schmidt(psi)
    #
    #evals, evecs = diag_lobpcg( pw, Vpot, psi, verbose=true, tol_avg=1e-10 )

    psi, Etot = Sch_solve_Emin_sd( pw, Vpot, psi, NiterMax=10 )
    psi, Etot = Sch_solve_Emin_cg( pw, Vpot, psi, NiterMax=1000 )
    #
    Y = ortho_gram_schmidt(psi)
    mu = Y' * op_H( pw, Vpot, Y )
    evals, evecs = eigen(mu)
    psi = Y*evecs

    for ist = 1:Nstates
        @printf("State # %d, Energy = %18.10f\n", ist, real(evals[ist]))
    end

    @printf("E_nn    = %18.10f\n", E_nn)
    @printf("E total = %18.10f\n", E_nn + Etot)

end

@time test_main( 64,64,64 )
