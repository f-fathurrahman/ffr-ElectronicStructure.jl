using Printf
using LinearAlgebra
using Random

include("PWGrid_v04.jl")
include("read_kpts.jl")
include("../common/gen_lattice.jl")

function test_main()
    LatConst = 10.
    LatVecs = gen_lattice_hexagonal( LatConst )
    ecutwfc_Ry = 20.0
    kpts_red = read_kpts("KPATH_HCP_60")

    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs, kpts_red )

    @printf("Number of G-vectors: %d\n", pw.gvec.Ng)
    @printf("Number of Gk-vectors\n")
    Ngw = pw.gkvec.Ngw
    for ik = 1:length(Ngw)
        @printf("%3d %8d\n", ik, Ngw[ik])
    end

    Ns = pw.Ns
    Npoints = prod(Ns) # note that for the current implementation Ns = Ng
    Ngwx = pw.gkvec.Ngwx

    @printf("Ngwx = %d\n", pw.gkvec.Ngwx)
    @printf("Ns   = (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])

end

test_main()
