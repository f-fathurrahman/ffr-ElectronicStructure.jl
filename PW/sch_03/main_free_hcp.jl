include("../pwgrid_04/PWGrid_v04.jl")
include("../pwgrid_04/read_kpts.jl")
include("../common/gen_lattice_ase.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("apply_K.jl")
include("apply_Vpot.jl")
include("apply_H.jl")
include("diag_lobpcg.jl")
include("Kprec.jl")

function test_main()
  const LatConst = 5.
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
  Nkpts = size(kpts_red)[2]
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
  for ik = 1:Nkpts
    @printf("%d %8.5f %8.5f %8.5f\n", ik, kpts_red[1,ik], kpts_red[2,ik], kpts_red[3,ik])
  end

  @printf("k-point list:\n")
  for ik = 1:Nkpts
    @printf("%d %8.5f %8.5f %8.5f\n", ik, kpts[1,ik], kpts[2,ik], kpts[3,ik])
  end

  #@printf("Number of Gk-vectors\n")
  #for ik = 1:length(Ngw)
  #  @printf("%3d %8d\n", ik, Ngw[ik])
  #end

  Vpot = zeros(Npoints)  # no external potential

  const Nstates = 5
  srand(1234)

  evals = zeros(Float64, Nstates, Nkpts )  # what is the optimal shape?

  for ik = 1:Nkpts

    Ngwk = Ngw[ik]
    psi  = randn(Ngwk,Nstates) + im*randn(Ngwk,Nstates)
    psi = ortho_gram_schmidt(psi)

    #Kpsi_k = apply_K( pw, psi, ik )
    #ss = sum(Kpsi_k)
    #@printf("Kin: %5d (%18.10f,%18.10f)\n", ik, real(ss), imag(ss))

    #Vpsi_k = apply_Vpot( pw, Vpot, psi, ik )
    #ss = sum(Vpsi_k)
    #@printf("Pot: %5d (%18.10f,%18.10f)\n", ik, real(ss), imag(ss))

    evals[:,ik], psi = diag_lobpcg( pw, Vpot, psi, ik, verbose=true, tol_avg=1e-7 )
  end

  plot_band_structure( evals, kpts )

end

using PyPlot
const plt = PyPlot

function plot_band_structure( evals, kpath )

  Nkpts = size(kpath)[2]
  Nstates = size(evals)[1]

  Xcoords = zeros(Float64,Nkpts)
  Xcoords[1] = 0.0
  for ik = 2:Nkpts
    dk = kpath[:,ik] - kpath[:,ik-1]
    Xcoords[ik] = Xcoords[ik-1] + norm( dk )
  end

  println(size(evals[1,:]'))
  println(size(Xcoords[:]))
  plt.clf()
  for is = 1:Nstates
    y = evals[is,:]'  # need transpose, UGH !
    plt.plot( Xcoords[:], y, marker="o")
  end
  plt.savefig("band_hcp_free.pdf")
end

test_main()