include("PsPot_HGH.jl")
include("ylmr2.jl")
include("../common/PWGrid_v02.jl")

function test_main()

  psp = PsPot_HGH(1, "Ge", "LDA_HGH/32ge.4.hgh")

  # for multispecies system this should be max over all species
  lmaxkb = psp.lmax
  @printf("lmaxkb = %d\n", lmaxkb)

  Ns = [32, 32, 32]
  LatVecs = 18.0*diagm(ones(3))
  pw = PWGrid( Ns, LatVecs )

  # use G-vectors for wavefunction
  G = pw.gvec.G[:,pw.gvecw.idx_gw2r]

  ylm = ylmr2( lmaxkb, G )

end

test_main()
