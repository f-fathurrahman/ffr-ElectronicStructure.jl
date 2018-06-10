using Printf
using LinearAlgebra

include("../common/PWGrid_v01.jl")
include("../common/gen_lattice.jl")
include("../common/print_matrix.jl")

function write_XSF( filnam, LL_in, atpos; molecule=false )
    #
    f = open(filnam, "w")
    Natoms = size(atpos)[2]
    #
    if molecule
        @printf(f, "MOLECULE\n")
    else
        @printf(f, "CRYSTAL\n")
    end
    LL = LL_in'  # Xcrysden uses lattice vectors by rows
    @printf(f, "PRIMVEC\n")
    @printf(f, "%18.10f %18.10f %18.10f\n", LL[1,1], LL[1,2], LL[1,3])
    @printf(f, "%18.10f %18.10f %18.10f\n", LL[2,1], LL[2,2], LL[2,3])
    @printf(f, "%18.10f %18.10f %18.10f\n", LL[3,1], LL[3,2], LL[3,3])
    @printf(f, "PRIMCOORD\n")
    @printf(f, "%8d %8d\n", Natoms, 1)
    for ia = 1:Natoms
        @printf(f, "X %18.10f %18.10f %18.10f\n", atpos[1,ia], atpos[2,ia], atpos[3,ia])
    end
    close(f)
end


function test_main_hexagonal()

    #LL = 16.0*diagm([1.0, 1.0, 1.0])
    #LL = gen_lattice_fcc(16.0)
    #LL = gen_lattice_bcc(16.0)

    Ns = [10, 10, 20]
    LL = gen_lattice_hexagonal(10.0, coa=2.0)
    pw = PWGrid( Ns, LL )
    atpos = pw.r

    println("Lattice vectors (by column)")
    print_matrix(pw.LatVecs)
    println("Reciprocal lattice vectors (by column)")    
    print_matrix(pw.RecVecs)

    write_XSF("R_grid_hexagonal.xsf", LL, atpos)

    Rec = pw.RecVecs*Ns[1]/2.0
    atpos = pw.G
    write_XSF("G_grid_hexagonal.xsf", LL, atpos, molecule=true)

    for ii = 1:3
        @printf("LatVecLen %d %18.10f\n", ii, norm(pw.LatVecs[:,ii]))
    end
    @printf("Ratio coa: %18.10f\n", norm(pw.LatVecs[:,1])/norm(pw.LatVecs[:,3]))

    @printf("\n")
    for ii = 1:3
        @printf("RecVecLen %d %18.10f\n", ii, norm(pw.RecVecs[ii,:]))
    end
    @printf("Ratio coa: %18.10f\n", norm(pw.RecVecs[:,3])/norm(pw.RecVecs[:,3]))

end


test_main_hexagonal()
