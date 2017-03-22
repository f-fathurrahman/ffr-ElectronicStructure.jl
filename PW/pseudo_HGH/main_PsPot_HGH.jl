include("PrintMatrix.jl")
include("PsPot_HGH.jl")

include("../common/PWGrid_v02.jl")

import PyPlot
const plt = PyPlot

function test_main()
  #psp1 = PsPot_HGH(1, "H" , "LDA_HGH/1h.1.hgh")
  #psp2 = PsPot_HGH(2, "Ti", "LDA_HGH/22ti.12.hgh")
  psp = PsPot_HGH(3, "Sm", "LDA_HGH/62sm.16.hgh")
  #psp4 = PsPot_HGH(4, "Pt", "LDA_HGH/78pt.18.hgh")

  info_PsPot_HGH(psp)

  Ns = [32, 32, 32]
  LatVecs = 5.0*diagm(ones(3))
  pw = PWGrid( Ns, LatVecs )

  # use G-vectors for wavefunction
  Gm = sqrt(pw.gvec.G2[pw.gvecw.idx_gw2r])
  Vprj_s = eval_HGH_proj_G( psp, 0, 1, Gm, pw.Ω )
  Vprj_p = eval_HGH_proj_G( psp, 1, 1, Gm, pw.Ω )
  Vprj_d = eval_HGH_proj_G( psp, 2, 1, Gm, pw.Ω )
  Vprj_f = eval_HGH_proj_G( psp, 3, 1, Gm, pw.Ω )

  plt.clf()
  plt.plot( Gm, Vprj_s, marker="o", label="s-prj" )
  plt.plot( Gm, Vprj_p, marker="o", label="p-prj" )
  plt.plot( Gm, Vprj_d, marker="o", label="d-prj" )
  plt.plot( Gm, Vprj_f, marker="o", label="f-prj" )
  plt.grid()
  plt.legend()
  plt.savefig("Sm_prj.png", dpi=300)


end

test_main()
