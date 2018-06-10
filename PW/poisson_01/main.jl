using Printf
using LinearAlgebra

include("../common/PWGrid_v01.jl")
include("../common/wrappers_fft.jl")
include("../common/gen_lattice.jl")

include("gen_dr.jl")
include("gen_rho.jl")
include("Poisson_solve.jl")

function test_main()
    #
    Ns = [64, 64, 64]
    LatVecs = gen_lattice_sc(16.0)
    #
    pw = PWGrid( Ns, LatVecs )
    #
    Npoints = pw.Npoints
    Ω = pw.Ω
    r = pw.r
    Ns = pw.Ns
    #
    # Generate array of distances
    #
    center = sum(LatVecs,dims=2)/2
    println("center = ", center)
    dr = gen_dr( r, center )
    #
    # Generate charge density
    #
    σ1 = 0.75
    σ2 = 0.50
    rho = gen_rho( dr, σ1, σ2 )
    #
    # Solve Poisson equation and calculate Hartree energy
    #
    phi = Poisson_solve( pw, rho )
    Ehartree = 0.5*dot( phi, rho ) * Ω/Npoints
    #
    Uanal = ( (1/σ1 + 1/σ2)/2 - sqrt(2) / sqrt( σ1^2 + σ2^2 ) ) / sqrt(pi)
    @printf("Num, ana, diff = %18.10f %18.10f %18.10e\n", Ehartree, Uanal, abs(Ehartree-Uanal))
end

@time test_main()
