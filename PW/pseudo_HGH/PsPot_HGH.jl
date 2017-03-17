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
function PsPot_HGH( itype, atsymb, filename )

  # Initialze with default values
  psp = PsPot_HGH( itype, atsymb, 0.0, 0, 0, 0.0, zeros(Float64,4), zeros(Float64,4),
                   zeros(Float64, 4,3,3), zeros(Float64, 4,3,3) )

  @printf("\nFile = %s\n", filename)
  file = open( filename )

  comment = readline(file)
  #@printf("%s", comment)

  lines = split( readline(file) )
  psp.zval = parse( Float64, lines[1] )
  zion = parse( Float64, lines[2] )
  #@printf("zion = %f\n", zion )
  #@printf("zval = %f\n", psp.zval )
  #@printf("pspdat = %s\n", lines[3] )

  lines  = split( readline(file) )
  pspcod = parse( Int, lines[1] )
  pspxc  = parse( Int, lines[2] )
  psp.lmax = parse( Int, lines[3] )
  psp.lloc = parse( Int, lines[4] )
  mmax   = parse( Int, lines[5] )
  r2well = parse( Int, lines[6] )
  #@printf("pspcod = %d\n", pspcod)
  #@printf("pspxc  = %d\n", pspxc)
  #@printf("lmax   = %d\n", psp.lmax)
  #@printf("lloc   = %d\n", psp.lloc)
  #@printf("mmax   = %d\n", mmax)
  #@printf("r2well = %d\n", r2well)

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
  #@printf("Angular momentum: %s\n", ANGMOM[l+1])
  #@printf("r: %f\n", psp.rc[l+1])
  #
  psp.h[1,1,1] = parse( lines[2] )
  psp.h[1,2,2] = parse( lines[3] )
  psp.h[1,3,3] = parse( lines[4] )
  #@printf("h: %f %f %f\n", psp.h[1,1,1], psp.h[1,2,2], psp.h[1,3,3])

  for l = 1:3
    lines = split( readline(file) )
    #
    psp.rc[l+1] = parse( lines[1] )
    #@printf("Angular momentum: %s\n", ANGMOM[l+1])
    #@printf("rc: %f\n", psp.rc[l+1])
    #
    psp.h[l+1,1,1] = parse( lines[2] )
    psp.h[l+1,2,2] = parse( lines[3] )
    psp.h[l+1,3,3] = parse( lines[4] )
    #@printf("h: %f %f %f\n", psp.h[l+1,1,1], psp.h[l+1,2,2], psp.h[l+1,3,3])
    #
    lines = split( readline(file) )
    psp.k[l+1,1,1] = parse( lines[1] )
    psp.k[l+1,2,2] = parse( lines[2] )
    psp.k[l+1,3,3] = parse( lines[3] )
    #
    #@printf("k: %f %f %f\n", psp.k[l+1,1,1], psp.k[l+1,2,2], psp.k[l+1,3,3])
  end

  close(file)

  const M_HALF = 0.5
  const M_THREE = 3.0
  const M_FIVE = 5.0
  const M_ONE = 1.0

  # from Octopus code
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
