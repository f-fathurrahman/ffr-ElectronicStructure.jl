using Printf
using LinearAlgebra
using Random

include("../pwgrid_04/PWGrid_v04.jl")
include("../pwgrid_04/read_kpts.jl")
include("../common/gen_lattice_ase.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("op_K.jl")
include("op_Vpot.jl")
include("op_H.jl")
include("diag_lobpcg.jl")
include("Kprec.jl")

include("structure_factor.jl")

include("dump_bandstructure.jl")

function test_main()
    LatConst = 5.0
    LatVecs = gen_lattice_hexagonal( LatConst )
    ecutwfc_Ry = 20.0
    kpts_red = read_kpts("../pwgrid_04/KPATH_HCP_60")

    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs, kpts_red )

    Ns = pw.Ns
    Npoints = prod(Ns) # note that for the current implementation Ns = Ng
    Ngwx = pw.gkvec.Ngwx
    LatVecs = pw.LatVecs
    RecVecs = pw.RecVecs
    Ngw = pw.gkvec.Ngw
    Nkpt = size(kpts_red)[2]
    kpts = pw.gkvec.kpts

    @printf("Unit cell:\n")
    for i = 1:3
        for j = 1:3
            @printf("%10.5f ", LatVecs[i,j])
        end
        @printf("\n")
    end

    @printf("Reciprocal unit cell:\n")
    for i = 1:3
        for j = 1:3
            @printf("%10.5f ", RecVecs[i,j])
        end
        @printf("\n")
    end

    @printf("Number of G-vectors: %d\n", pw.gvec.Ng)
    @printf("Ngwx = %d\n", pw.gkvec.Ngwx)
    @printf("Ns   = (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])

    @printf("k-point list: (reduced)\n")
    for ik = 1:Nkpt
        @printf("%d %8.5f %8.5f %8.5f\n", ik, kpts_red[1,ik], kpts_red[2,ik], kpts_red[3,ik])
    end

    @printf("k-point list:\n")
    for ik = 1:Nkpt
        @printf("%d %8.5f %8.5f %8.5f\n", ik, kpts[1,ik], kpts[2,ik], kpts[3,ik])
    end

    Xpos = reshape( [0.0, 0.0, 0.0], (3,1) )
    Sf = structure_factor( Xpos, pw.gvec.G )
    Vg = zeros(ComplexF64,Npoints)
    prefactor = -4*pi/pw.Ω
    for ig=2:Npoints
        Vg[ig] = prefactor/pw.gvec.G2[ig]
    end
    Vpot = real( G_to_R(Ns, Vg .* Sf) ) * Npoints

    Nstates = 4
    srand(1234)

    evals = zeros(Float64, Nstates, Nkpt)  # what is the optimal shape?

    for ik = 1:Nkpt

        Ngwk = Ngw[ik]
        psi  = rand(ComplexF64,Ngwk,Nstates)
        psi = ortho_gram_schmidt(psi)

        evals[:,ik], psi = diag_lobpcg( pw, Vpot, psi, ik, verbose=true, tol_avg=1e-7 )
    end

    dump_bandstructure( evals, kpts, filename="band_H_hcp.dat" )

end


test_main()
