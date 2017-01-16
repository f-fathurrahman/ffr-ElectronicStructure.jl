type PsPot_HGH
  itype::Int
  atsymb::ASCIIString
  zval::Float64
  lloc::Int
  lmax::Int
  rloc::Float64
  rc::Array{Float64}   # indexed (l+1), l=0,1,2,3
  c::Array{Float64}    # indexed 1,2,3,4
  h::Array{Float64,3}  # indexed [0:3,1:3,1:3]
  k::Array{Float64,3}  # indexed [0:3,1:3,1:3]
end


function PsPot_HGH( itype, atsymb, filename )

  # Default values
  psp = PsPot_HGH( itype, atsymb, 0.0, 0, 0, 0.0, zeros(Float64,4), zeros(Float64,4),
                   zeros(Float64, 4,3,3), zeros(Float64, 4,3,3) )

  @printf("\nFile = %s\n", filename)
  file = open( filename )

  comment = readline(file)
  #@printf("%s", comment)

  lines = split( readline(file) )
  psp.zval = parse( Float64, lines[1] )
  zion = parse( Float64, lines[2] )
  @printf("zion = %f\n", zion )
  @printf("zval = %f\n", psp.zval )
  @printf("pspdat = %s\n", lines[3] )

  lines = split( readline(file) )
  pspcod = parse( Int, lines[1] )
  pspxc  = parse( Int, lines[2] )
  psp.lmax   = parse( Int, lines[3] )
  psp.lloc   = parse( Int, lines[4] )
  mmax   = parse( Int, lines[5] )
  r2well = parse( Int, lines[6] )
  @printf("pspcod = %d\n", pspcod)
  @printf("pspxc  = %d\n", pspxc)
  @printf("lmax   = %d\n", psp.lmax)
  @printf("lloc   = %d\n", psp.lloc)
  @printf("mmax   = %d\n", mmax)
  @printf("r2well = %d\n", r2well)

  lines = split( readline(file) )
  rloc = parse( Float64, lines[1] )
  psp.c[1] = parse( Float64, lines[2] )
  psp.c[2] = parse( Float64, lines[3] )
  psp.c[3] = parse( Float64, lines[4] )
  psp.c[4] = parse( Float64, lines[5] )

  const ANGMOM = ["s", "p", "d", "f"]

  l = 0   # s
  lines = split( readline(file) )
  #
  psp.rc[l+1] = parse( lines[1] )
  @printf("Angular momentum: %s\n", ANGMOM[l+1])
  @printf("r: %f\n", psp.rc[l+1])
  #
  psp.h[l+1,1,1] = parse( lines[1] )
  psp.h[l+1,2,2] = parse( lines[2] )
  psp.h[l+1,2,2] = parse( lines[3] )
  @printf("h: %f %f %f\n", psp.h[l+1,1,1], psp.h[l+1,2,2], psp.h[l+1,3,3])

  for l = 1:3
    lines = split( readline(file) )
    #
    psp.rc[l+1] = parse( lines[1] )
    @printf("Angular momentum: %s\n", ANGMOM[l+1])
    @printf("rc: %f\n", psp.rc[l+1])
    #
    psp.h[l+1,1,1] = parse( lines[1] )
    psp.h[l+1,2,2] = parse( lines[2] )
    psp.h[l+1,2,2] = parse( lines[3] )
    @printf("h: %f %f %f\n", psp.h[l+1,1,1], psp.h[l+1,2,2], psp.h[l+1,3,3])
    #
    lines = split( readline(file) )
    psp.k[l+1,1,1] = parse( lines[1] )
    psp.k[l+1,2,2] = parse( lines[2] )
    psp.k[l+1,2,2] = parse( lines[3] )
    #
    @printf("k: %f %f %f\n", psp.k[l+1,1,1], psp.k[l+1,2,2], psp.k[l+1,3,3])
  end

  close(file)

  return psp
end


function test_main()
  psp1 = PsPot_HGH(1, "H" , "LDA_HGH/1h.1.hgh")
  psp2 = PsPot_HGH(2, "Ti", "LDA_HGH/22ti.12.hgh")
  psp3 = PsPot_HGH(3, "Sm", "LDA_HGH/62sm.16.hgh")
end

test_main()
