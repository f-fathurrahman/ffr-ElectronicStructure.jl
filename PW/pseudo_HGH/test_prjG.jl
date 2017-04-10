include("PsPot_HGH.jl")
include("ylmr2.jl")
include("../common/PWGrid_v02.jl")

function test_main()

  const Nspecies = 2  # test for one species only
  psp = Array{PsPot_HGH,1}(Nspecies)

  #psp[1] = PsPot_HGH(1, "Ge", "LDA_HGH/32ge.4.hgh")

  psp[1] = PsPot_HGH(1, "Cu", "LDA_HGH/29cu.11.hgh")
  psp[2] = PsPot_HGH(2, "Mo", "LDA_HGH/42mo.14.hgh")

  # for multispecies system this should be max over all species
  lmaxkb = 0
  for isp = 1:Nspecies
    if lmaxkb < psp[isp].lmax
      lmaxkb = psp[isp].lmax
    end
  end
  @printf("lmaxkb = %d\n", lmaxkb)

  Ns = [32, 32, 32]
  LatVecs = 18.0*diagm(ones(3))
  pw = PWGrid( Ns, LatVecs )
  Ω = pw.Ω

  # use G-vectors for wavefunction
  G = pw.gvec.G[:,pw.gvecw.idx_gw2r]
  Gm = sqrt(pw.gvec.G2[pw.gvecw.idx_gw2r])

  Ngwx = pw.gvecw.Ngwx

  ylm = ylmr2( lmaxkb, G )

  # from init_us_1

  NH = zeros(Int,Nspecies)

  for isp = 1:Nspecies
    ih = 0
    println("snprj = ", psp[isp].snprj)
    for ib = 1:psp[isp].snprj
      l = psp[isp].lll[ib]
      println("l = ", l)
      for m = 1:(2*l+1)
        #nhtol[ih,isp] = l
        #nhtolm[ih,isp] = l*l + m
        #indv[ih,isp] = ib
        ih = ih + 1
      end
    end
    NH[isp] = ih
    @printf("isp = %d, ih = %d\n", isp, ih)
  end

  Nhmax = maximum(NH)
  Nkb = sum(NH)
  println("Nhmax = ", Nhmax)
  println("Nkb   = ", Nkb)

  nhtol  = zeros(Int,Nhmax,Nspecies)
  nhtolm = zeros(Int,Nhmax,Nspecies)
  indv   = zeros(Int,Nhmax,Nspecies)

  for isp = 1:Nspecies
    ih = 1
    for ib = 1:psp[isp].snprj
      l = psp[isp].lll[ib]
      for m = 1:(2*l+1)
        nhtol[ih,isp] = l
        nhtolm[ih,isp] = l*l + m
        indv[ih,isp] = ib
        ih = ih + 1
      end
    end
  end

  vkb1 = zeros(Float64,Ngwx,Nhmax)
  Vkb = zeros(Complex128,Ngwx,Nkb)

  for isp = 1:Nspecies
    snprj = psp[isp].snprj
    @printf("\nSpecies = %d\n", isp)
    for ib = 1:snprj
      l = psp[isp].lll[ib]
      iproj = psp[isp].ipr[ib]
      # remember that Gm (the magnitudes) are used, not G
      Vprj = eval_HGH_proj_G( psp[isp], l, iproj, Gm, Ω )
      @printf("l = %d, iproj = %d done\n", l, iproj)

      for ih = 1:NH[isp]
        if ib == indv[ih,isp]
          lm = nhtolm[ih,isp]
          @printf("Same ib = %d, lm = %d\n", ib, lm)
          for ig=1:Ngwx
            vkb1[ig,ih] = ylm[ig,lm]*Vprj[ig]
          end
        end
      end

    end
  end

end

test_main()
