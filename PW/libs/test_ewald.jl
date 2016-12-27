include("../pwgrid_02/PWGrid.jl")

function latvec_hexagonal(a; ca=1.0)
  LL = zeros(3,3)
  LL[1,:] = [1.0, 0.0, 0.0]
  LL[2,:] = [cos(pi/3.0), sin(pi/3.0), 0.0]
  LL[3,:] = [0.0, 0.0, ca]
  return a*LL
end

function ewald()

  Ns = [60, 60, 60]
  const alat = 16.0
  LatVecs = alat*diamg(ones(3))
  pw = PWGrid( Ns, LatVecs )

  atpos = [0.0, 0.0, 0.0]

  nat = 1
  ntyp = 1
  ityp = [1]  # ityp(nat)
  zv = [1.0]     # zv(ntyp)  # need to conform with `ntypx` in PWSCF ?
  at = LatVecs'/alat
  bg = inv(at')
  #bg = pw.RecVecs'/(2*pi) # ???
  tau = atpos/alat
  omega = pw.Ω
  g  = pw.G
  gg = pw.G2
  ngm = prod(Ns) # for the moment
  gcutm =

  println( at )
  println( bg )
  println( at * bg' )

  type_args = (
               Ref{Float64},    # alat
               Ref{Int64},      # nat
               Ref{Int64},      # ntyp
               Ptr{Int64},      # ityp
               Ptr{Float64},    # zv
               Ptr{Float64},    # at
               Ptr{Float64},    # bg
               Ptr{Float64},    # tau
               Ref{Float64},    # omega
               Ptr{Float64},    # g
               Ptr{Float64},    # gg
               Ref{Int64},      # ngm
               Ref{Float64},    # gcutm
               Ref{Int64},      # gstart
               Ref{Int32},      # gamma_only
               Ptr{Complex128}, # strf
               )

  #ccall( (:ewald_, "libqe.so"), Float64,
  #       type_args,
  #       alat, nat, ntyp, ityp, zv, at, bg, tau, omega, g, gg,
  #       ngm, gcutm, gstart, gamma_only, strf )
end


ewald()
