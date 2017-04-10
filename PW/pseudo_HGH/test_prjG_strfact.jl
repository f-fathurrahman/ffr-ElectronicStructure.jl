include("PsPot_HGH.jl")
include("ylmr2.jl")
include("../common/PWGrid_v02.jl")
include("../common/calc_strfact_v2.jl")

function test_main()

  #
  # Setup plane wave basis
  #
  Ns = [32, 32, 32]
  LatVecs = 18.0*diagm(ones(3))
  pw = PWGrid( Ns, LatVecs )
  Ω = pw.Ω

  # use G-vectors for wavefunction
  # CAUTION: this should be depend on k-points
  idx_gw2r = pw.gvecw.idx_gw2r
  G = pw.gvec.G[:,idx_gw2r]
  Gm = sqrt( pw.gvec.G2[idx_gw2r] )

  Ngwx = pw.gvecw.Ngwx
  Ng = pw.gvec.Ng

  #
  # Setup structure and pseudopotentials
  #
  const Nspecies = 2
  const Natoms   = 3
  atmsymb = ["Cu", "Mo"] # unique list of atomic symbols
  atm2species = [1, 2, 1]  # mapping from atom to species

  # not really used
  Zv = zeros(Float64,Nspecies)
  Zv = [11.0, 14.0]

  Xpos = reshape( [0.0, 0.0, 0.0,  # Cu
                   1.5, 0.0, 0.0,  # Mo
                   0.0, 1.5, 0.0], # Cu
                   (3,Natoms) )

  # Calculate structure factor
  Sf = calc_strfact( Xpos, Nspecies, atm2species, pw.gvec.G ) # Sf[Ng,Nspecies]

  psp = Array{PsPot_HGH,1}(Nspecies)

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

  println("Nhmax = ", Nhmax)

  Nkb = 0
  for ia = 1:Natoms
    isp = atm2species[ia]
    Nkb = Nkb + NH[isp]
  end
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

  jkb = 0

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
          # XXXXXX  need to zero out vkb1 ???
        end
      end
    end

    for ia = 1:Natoms
      if atm2species[ia] == isp
        # arg = dot( kpts[:], atpos[:,ia] ) )  # needed for kpts other than gamma (0,0,0)
        phase = 1.   # special case for gamma point sampling only
        for ih = 1:NH[isp]
          jkb = jkb + 1
          pre = (-1.0*im)^nhtol[ih,isp]*phase
          for ig = 1:Ngwx
            Vkb[ig,jkb] = vkb1[ig,ih]*Sf[idx_gw2r[ig],isp]*pre
          end
          # need to zero out multiple k-points
        end
      end
    end

  end # loop over species

end

test_main()
