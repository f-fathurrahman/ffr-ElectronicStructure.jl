using Printf
using LinearAlgebra

include("../common/PWGrid_v03.jl")
include("../common/gen_lattice_pwscf.jl")

function test_main()
    LatVecs = gen_lattice_sc(8.0)

    pw = PWGrid( 15.0, LatVecs )
    Ns = pw.Ns
    @printf("Sampling points: (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])

    actual = prod(pw.Ns)/pw.gvecw.Ngwx
    theor = 1/(4*pi*0.25^3/3)
    @printf("Actual, theor: %10.5f %10.5f\n", actual, theor)

end

@time test_main()
