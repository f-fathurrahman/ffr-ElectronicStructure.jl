using Printf
using LinearAlgebra
using Random

include("../common/PWGrid_v01.jl")
include("../common/wrappers_fft.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/gen_lattice_pwscf.jl")

include("gen_dr.jl")
include("init_pot_harm_3d.jl")
include("op_K.jl")
include("op_Vpot.jl")
include("op_H.jl")

include("calc_rho.jl")
include("calc_Etot.jl")
include("Kprec.jl")
include("diag_lobpcg.jl")

function test_main( ns1::Int,ns2::Int,ns3::Int )
    #
    Ns = [ns1,ns2,ns3]
    LatVecs = gen_lattice_sc(6.0)
    #
    pw = PWGrid( Ns, LatVecs )
    #
    Npoints = pw.Npoints
    Ω = pw.Ω
    r = pw.r
    G = pw.G
    G2 = pw.G2
    #
    # Generate array of distances
    #
    center = sum(LatVecs,dims=2)/2
    dr = gen_dr( r, center )
    #
    # Setup potential
    #
    Vpot = init_pot_harm_3d( pw, dr )
    println("sum(Vpot)*Ω/Npoints = ", sum(Vpot)*Ω/Npoints);
    #
    Nstates = 4
    Random.seed!(2222)
    psi = randn(ComplexF64,Npoints,Nstates)
    psi = ortho_gram_schmidt(psi)
    #
    evals, evecs = diag_lobpcg( pw, Vpot, psi, verbose=true, tol_avg=1e-10 )

    for ist = 1:Nstates
        @printf("State %3d, Energy = %18.10f\n", ist, real(evals[ist]))
    end
end

@time test_main( 30,30,30 )
