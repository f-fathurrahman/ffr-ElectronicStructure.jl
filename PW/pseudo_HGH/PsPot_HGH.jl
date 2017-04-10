type PsPot_HGH
  itype::Int
  atsymb::ASCIIString
  zval::Float64
  lloc::Int
  lmax::Int
  rloc::Float64
  rc::Array{Float64}   # indexed (l+1), l=0,1,2,3
  c::Array{Float64}    # indexed 1,2,3,4
  h::Array{Float64,3}  # originally indexed [0:3,1:3,1:3]
  k::Array{Float64,3}  # indexed [0:3,1:3,1:3]
end


# Constructor
function PsPot_HGH( itype::Int, atsymb::ASCIIString, filename::ASCIIString )

  # Initialze with default values
  psp = PsPot_HGH( itype, atsymb, 0.0, 0, 0, 0.0, zeros(Float64,4), zeros(Float64,4),
                   zeros(Float64, 4,3,3), zeros(Float64, 4,3,3) )

  @printf("\nFile = %s\n", filename)
  file = open( filename )

  comment = readline(file)

  lines = split( readline(file) )
  psp.zval = parse( Float64, lines[1] )
  zion = parse( Float64, lines[2] )

  lines  = split( readline(file) )
  pspcod = parse( Int, lines[1] )
  pspxc  = parse( Int, lines[2] )
  psp.lmax = parse( Int, lines[3] )
  psp.lloc = parse( Int, lines[4] )
  mmax   = parse( Int, lines[5] )
  r2well = parse( Int, lines[6] )

  lines = split( readline(file) )
  psp.rloc = parse( Float64, lines[1] )
  psp.c[1] = parse( Float64, lines[2] )
  psp.c[2] = parse( Float64, lines[3] )
  psp.c[3] = parse( Float64, lines[4] )
  psp.c[4] = parse( Float64, lines[5] )

  const ANGMOM = ["s", "p", "d", "f"]

  l = 0  # s
  lines = split( readline(file) )
  #
  psp.rc[1] = parse( lines[1] )
  #
  psp.h[1,1,1] = parse( lines[2] )
  psp.h[1,2,2] = parse( lines[3] )
  psp.h[1,3,3] = parse( lines[4] )

  for l = 1:3
    lines = split( readline(file) )
    #
    psp.rc[l+1] = parse( lines[1] )
    #
    psp.h[l+1,1,1] = parse( lines[2] )
    psp.h[l+1,2,2] = parse( lines[3] )
    psp.h[l+1,3,3] = parse( lines[4] )
    #
    lines = split( readline(file) )
    psp.k[l+1,1,1] = parse( lines[1] )
    psp.k[l+1,2,2] = parse( lines[2] )
    psp.k[l+1,3,3] = parse( lines[3] )
    #
  end

  close(file)

  const M_HALF = 0.5
  const M_THREE = 3.0
  const M_FIVE = 5.0
  const M_ONE = 1.0

  # from Octopus code, see also the appendix of HGH paper

  psp.h[0+1, 1, 2] = -M_HALF    * sqrt(M_THREE/M_FIVE) * psp.h[0+1, 2, 2]
  psp.h[0+1, 1, 3] =  M_HALF    * sqrt(M_FIVE/21.0)    * psp.h[0+1, 3, 3]
  psp.h[0+1, 2, 3] = -M_HALF    * sqrt(100.0/63.0)     * psp.h[0+1, 3, 3]
  psp.h[1+1, 1, 2] = -M_HALF    * sqrt(M_FIVE/7.0)     * psp.h[1+1, 2, 2]
  psp.h[1+1, 1, 3] =  M_ONE/6.0 * sqrt(35.0/11.0)      * psp.h[1+1, 3, 3]
  psp.h[1+1, 2, 3] = -M_ONE/6.0 * ( 14.0 / sqrt(11.0)) * psp.h[1+1, 3, 3]
  psp.h[2+1, 1, 2] = -M_HALF    * sqrt(7.0/9.0)        * psp.h[2+1, 2, 2]
  psp.h[2+1, 1, 3] =  M_HALF    * sqrt(63.0/143.0)     * psp.h[2+1, 3, 3]
  psp.h[2+1, 2, 3] = -M_HALF    * (18.0/sqrt(143.0))   * psp.h[2+1, 3, 3]

  psp.k[0+1, 1, 2] = -M_HALF    * sqrt(M_THREE/M_FIVE) * psp.k[0+1, 2, 2]
  psp.k[0+1, 1, 3] =  M_HALF    * sqrt(M_FIVE/21.0)    * psp.k[0+1, 3, 3]
  psp.k[0+1, 2, 3] = -M_HALF    * sqrt(100.0/63.0)     * psp.k[0+1, 3, 3]
  psp.k[1+1, 1, 2] = -M_HALF    * sqrt(M_FIVE/7.0)     * psp.k[1+1, 2, 2]
  psp.k[1+1, 1, 3] =  M_ONE/6.0 * sqrt(35.0/11.0)      * psp.k[1+1, 3, 3]
  psp.k[1+1, 2, 3] = -M_ONE/6.0 * (14.0 / sqrt(11.0))  * psp.k[1+1, 3, 3]
  psp.k[2+1, 1, 2] = -M_HALF    * sqrt(7.0/9.0)        * psp.k[2+1, 2, 2]
  psp.k[2+1, 1, 3] =  M_HALF    * sqrt(63.0/143.0)     * psp.k[2+1, 3, 3]
  psp.k[2+1, 2, 3] = -M_HALF    * (18.0/sqrt(143.0))   * psp.k[2+1, 3, 3]

  # Parameters are symmetric.
  for k = 1:4
    for i = 1:3
      for j = i+1:3
        psp.h[k, j, i] = psp.h[k, i, j]
        psp.k[k, j, i] = psp.k[k, i, j]
      end
    end
  end

  return psp
end


function info_PsPot_HGH( psp::PsPot_HGH )

  const ANGMOM = ["s", "p", "d", "f"]

  @printf("rloc: %f, c: %f, %f, %f, %f\n", psp.rloc, psp.c[1], psp.c[2], psp.c[3], psp.c[4])
  for k=1:4
    @printf("Angular momentum: %s, rc = %f\n", ANGMOM[k], psp.rc[k])
    @printf("h = \n")
    PrintMatrix( reshape(psp.h[k,:,:],(3,3) ) )
  end

end

# Evaluate HGH projector function in G-space
function eval_HGH_Vloc_G( psp, G2, Ω )

  Ng = size(G2)[1]
  Vg = zeros(Ng)

  rloc = psp.rloc
  zval = psp.zval
  c1 = psp.c[1]
  c2 = psp.c[2]
  c3 = psp.c[3]
  c4 = psp.c[4]

  pre1 = -4*pi*zval/Ω
  pre2 = sqrt(8*pi^3)*rloc^3/Ω
  #
  for ig=2:Ng
    Gr = sqrt(G2[ig])*rloc
    expGr2 = exp(-0.5*Gr^2)
    Vg[ig] = pre1/G2[ig]*expGr2 + pre2*expGr2 * (c1 + c2*(3-Gr^2) +
             c3*(15 - 10*Gr^2 + Gr^4) + c4*(105 - 105*Gr^2 + 21*Gr^4 - Gr^6) )
  end
  # limiting value, with minus sign ?
  Vg[1] = 2*pi*zval*rloc^2 + (2*pi)^1.5 * rloc^3 * (c1 + 3.0*c2 + 15*c3 + 105*c4)

  return Vg
end


# Evaluate HGH projector function in G-space
function eval_HGH_proj_G( psp, l, iproj, G, Ω )

  # G is magnitude of G-vectors
  Ng = size(G)[1]

  Vprj = zeros(Ng)

  rrl = psp.rc[l+1]

  # s-channel
  if l == 0

    if iproj==1

      for ig = 1:Ng
        Gr2 = ( G[ig]*rrl )^2
        Vprj[ig] = exp( -0.5*Gr2 )
      end

    elseif iproj==2

      for ig = 1:Ng
        Gr2 = ( G[ig]*rrl )^2
        Vprj[ig] = 2.0/sqrt(15.0) * exp( -0.5*Gr2 ) * ( 3.0 - Gr2 )
      end

    elseif iproj==3

      for ig = 1:Ng
        Gr2 = ( G[ig]*rrl )^2
        Vprj[ig] = (4.0/3.0)/sqrt(105.0) * exp( -0.5*Gr2 ) * (15.0 - 10.*Gr2 + Gr2^2)
       end

    end  # if iproj

  # p-channel
  elseif l == 1

    if iproj == 1

      for ig = 1:Ng
        Gr2 = ( G[ig]*rrl )^2
        Vprj[ig] = (1.0/sqrt(3.0)) * exp(-0.5*Gr2) * G[ig]
      end

    elseif iproj == 2

      for ig = 1:Ng
        Gr2 = (G[ig]*rrl)^2
        Vprj[ig] = (2./sqrt(105.)) * exp(-0.5*Gr2) * G[ig]*(5. - Gr2)
      end

    elseif iproj == 3

      for ig = 1:Ng
        Gr2 = ( G[ig]*rrl)^2
        Vprj[ig] = (4./3.)/sqrt(1155.) * exp(-0.5*Gr2) * G[ig] * (35. - 14.*Gr2 + Gr2^2)
      end

    end # if iproj

  # d-channel
  elseif l == 2

    if iproj == 1

      for ig = 1:Ng
        Gr2 = ( G[ig]*rrl )^2
        Vprj[ig] = (1.0/sqrt(15.0)) * exp(-0.5*Gr2) * G[ig]^2
      end

    elseif iproj == 2

      for ig = 1:Ng
        Gr2 = (G[ig]*rrl)^2
        Vprj[ig] = (2./3.)/sqrt(105.) * exp(-0.5*Gr2) * G[ig]^2 * (7.-Gr2)
      end

    end # if iproj

  # f-channel
  elseif l == 3

    for ig = 1:Ng
      Gr2 = ( G[ig]*rrl )^2
      Vprj[ig] = G[ig]^3 * exp(-0.5*Gr2)
    end

  end  # if l

  pre =  4 * pi^(5./4.) * sqrt( 2.^(l+1) * rrl^(2*l+3) / Ω )

  Vprj[:] = pre * Vprj[:]

  return Vprj

end
