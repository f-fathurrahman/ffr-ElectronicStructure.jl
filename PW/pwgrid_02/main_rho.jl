include("../common/PWGrid_v01.jl")
include("../common/wrappers_fft_v01.jl")

include("structure_factor.jl")
include("gen_dr.jl")
include("gen_rho.jl")
include("write_xsf.jl")

function test_main()
    Ns = [64,64,64]
    const a = 16.0
    LatVecs = a*diagm( [1.0, 1.0, 1.0] )
    pwgrid = PWGrid(Ns,LatVecs)

    Npoints = prod(Ns)
    Ω = pwgrid.Ω
    G2 = pwgrid.G2

    # Atomic positions and nuclear charge
    #Xpos = reshape( [5.0, 5.0, 5.0,
    #                 6.5, 5.0, 5.0], (3,2) )
    Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )

    # Calculate structure factor
    @printf("Npoints = %8d\n", Npoints)
    Sf = structure_factor( Xpos, pwgrid.G )

    center = 0.5*a*ones(3)
    println(center)
    dr = gen_dr( pwgrid.r, center )
    rho = gen_rho( Ns, dr, 0.25, Sf )

    sumSf = sum(Sf)
    @printf("sum(Sf) = (%15.10e,%15.10e)\n", real(sumSf), imag(sumSf))
    @printf("integRho = %18.10f\n", sum(rho)*Ω/Npoints)

    write_xsf( "rho.xsf", LatVecs, Xpos )
    write_xsf_3d_crystal( "rho.xsf", Ns, LatVecs, rho )
end

@time test_main()
