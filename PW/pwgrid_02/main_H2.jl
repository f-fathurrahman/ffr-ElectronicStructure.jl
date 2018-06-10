using Printf
using LinearAlgebra

include("../common/PWGrid_v01.jl")
include("../common/wrappers_fft.jl")
include("../common/gen_lattice_pwscf.jl")
include("structure_factor.jl")
include("write_xsf.jl")

function test_main()
    Ns = [50,50,50]
    LatVecs = gen_lattice_sc(10.0)
    pwgrid = PWGrid(Ns,LatVecs)

    Npoints = prod(Ns)
    Ω = pwgrid.Ω
    G2 = pwgrid.G2

    # Atomic positions and nuclear charge
    Xpos = reshape( [5.0, 5.0, 5.0,
                     6.5, 5.0, 5.0], (3,2) )
    #Xpos = reshape( [4.0, 4.0, 4.0], (3,1) )

    # Calculate structure factor
    Sf = structure_factor( Xpos, pwgrid.G )

    Ng = Npoints
    Vg = zeros(ComplexF64, Ng)
    prefact = -4.0*pi/Ω
    for ig = 2:Ng
        Vg[ig] = prefact/G2[ig]
    end
    Vg = Vg .* Sf
    Vr = real( G_to_R(Ns, Vg) )*Npoints

    @printf("sum(Vr) = %18.10e\n", sum(Vr))

    write_xsf( "H2.xsf", LatVecs, Xpos )
    write_xsf_3d_crystal( "H2.xsf", Ns, LatVecs, Vr )
end

@time test_main()
