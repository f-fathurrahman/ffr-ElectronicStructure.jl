! NOTE: This will only works for orthorombic cell
SUBROUTINE calc_dr_periodic( LL, r0, Npoints, lingrid, dr )
  IMPLICIT NONE
  !
  REAL(8) :: LL(3)
  REAL(8) :: r0(3)
  INTEGER :: Npoints
  REAL(8) :: lingrid(3,Npoints)
  REAL(8) :: dr(Npoints)
  !
  INTEGER :: ip
  REAL(8) :: Lx, Ly, Lz
  REAL(8) :: xx1, xx2, xx3, xx
  REAL(8) :: yy1, yy2, yy3, yy
  REAL(8) :: zz1, zz2, zz3, zz

  Lx = LL(1)
  Ly = LL(2)
  Lz = LL(3)

  DO ip = 1, Npoints
    xx1 = abs( lingrid(1,ip) - r0(1) )
    xx2 = abs( r0(1) + Lx - lingrid(1,ip) )
    xx3 = abs( lingrid(1,ip) - r0(1) - Lx)
    xx  = minval( (/xx1,xx2,xx3/) )
    !
    yy1 = abs( lingrid(2,ip) - r0(2) )
    yy2 = abs( r0(2) + Ly - lingrid(2,ip) )
    yy3 = abs( lingrid(2,ip) - r0(2) - Ly)
    yy  = minval( (/yy1,yy2,yy3/) )
    !
    zz1 = abs( lingrid(3,ip) - r0(3) )
    zz2 = abs( r0(3) + Lz - lingrid(3,ip) )
    zz3 = abs( lingrid(3,ip) - r0(3) - Lz)
    zz  = minval( (/zz1,zz2,zz3/) )
    !
    dr(ip) = sqrt( xx**2 + yy**2 + zz**2 )
  ENDDO

END SUBROUTINE 
