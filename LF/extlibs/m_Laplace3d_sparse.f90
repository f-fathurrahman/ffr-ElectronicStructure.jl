! efefer 1 January 2016
! modified 20 October 2016

SUBROUTINE init_Laplacian3d( NNZ, NN, Npoints, lin2xyz, x_D2jl, y_D2jl, z_D2jl, row, column, values )
  IMPLICIT NONE
  !
  INTEGER :: NNZ
  INTEGER(8) :: NN(3)
  INTEGER(8) :: Npoints
  INTEGER(8) :: lin2xyz(3,Npoints)
  REAL(8) :: x_D2jl( NN(1), NN(1) )
  REAL(8) :: y_D2jl( NN(2), NN(2) )
  REAL(8) :: z_D2jl( NN(3), NN(3) )
  INTEGER(8) :: row(NNZ), column(NNZ)
  REAL(8) :: values(NNZ)
  !
  INTEGER(8) :: ip1, ip2, i1,i2, j1,j2, k1,k2
  REAL(8) :: nabla2
  LOGICAL :: Tnz(3)
  INTEGER :: inz

  !WRITE(*,*) 'Calling init_Laplacian'
  !WRITE(*,*) 'NN = ', NN(1:3)
  !WRITE(*,*) 'Npoints = ', Npoints

  !WRITE(*,*) 'Pass here 46'
  inz = 0
  !WRITE(*,*) 'Pass here 48'
  DO ip2 = 1, Npoints
    DO ip1 = ip2, Npoints
      Tnz(:) = .FALSE.
      i1 = lin2xyz(1,ip1)
      i2 = lin2xyz(1,ip2)
      !
      j1 = lin2xyz(2,ip1)
      j2 = lin2xyz(2,ip2)
      !
      k1 = lin2xyz(3,ip1)
      k2 = lin2xyz(3,ip2)
      !
      nabla2 = 0.d0
      !
      IF( j1 == j2 .AND. k1 == k2 ) THEN
        nabla2 = nabla2 + x_D2jl(i1,i2)
        Tnz(1) = .TRUE.
      ENDIF
      !
      IF( i1 == i2 .AND. k1 == k2 ) THEN
        nabla2 = nabla2 + y_D2jl(j1,j2)
        Tnz(2) = .TRUE.
      ENDIF
      !
      IF( i1 == i2 .AND. j1 == j2 ) THEN
        nabla2 = nabla2 + z_D2jl(k1,k2)
        Tnz(3) = .TRUE.
      ENDIF
      !WRITE(*,*) ip1,ip2, any(Tnz)
      IF( ANY(Tnz) ) THEN
        inz = inz + 1
        values(inz) = nabla2
        column(inz) = ip1
        row(inz)    = ip2
        !WRITE(*,*) 'Nonzero:', ip1, ip2, nabla2
        IF( ip1 /= ip2 ) THEN
          inz = inz + 1
          values(inz) = nabla2
          column(inz) = ip2
          row(inz)    = ip1
        ENDIF
      ENDIF ! Tnz
      !
    ENDDO !
  ENDDO

  !WRITE(*,*) 'NNZ, Npoints = ', NNZ, Npoints
  !WRITE(*,*) 'Percentage   = ', ( NNZ ) / dble(Npoints**2) * 100d0
END SUBROUTINE


FUNCTION get_NNZ_Laplacian3d( NN, Npoints, lin2xyz, x_D2jl, y_D2jl, z_D2jl ) RESULT(NNZ)
  IMPLICIT NONE
  !
  INTEGER(8) :: NN(3)
  INTEGER(8) :: Npoints
  INTEGER(8) :: lin2xyz(3,Npoints)
  REAL(8) :: x_D2jl( NN(1), NN(1) )
  REAL(8) :: y_D2jl( NN(2), NN(2) )
  REAL(8) :: z_D2jl( NN(3), NN(3) )
  INTEGER(8) :: NNZ
  !
  INTEGER(8) :: istep, ip1, ip2, i1,i2, j1,j2, k1,k2
  REAL(8) :: nabla2
  LOGICAL :: Tnz(3)

  !WRITE(*,*) 'Calling init_Laplacian'
  !WRITE(*,*) 'NN = ', NN(1:3)
  !WRITE(*,*) 'Npoints = ', Npoints

  !WRITE(*,*) 'Pass here 46'
  NNZ = 0
  !WRITE(*,*) 'Pass here 48'
  DO ip2 = 1, Npoints
    DO ip1 = ip2, Npoints
      !WRITE(*,*) ip1, ip2
      Tnz(:) = .FALSE.
      i1 = lin2xyz(1,ip1)
      i2 = lin2xyz(1,ip2)
      !
      j1 = lin2xyz(2,ip1)
      j2 = lin2xyz(2,ip2)
      !
      k1 = lin2xyz(3,ip1)
      k2 = lin2xyz(3,ip2)
      !
      nabla2 = 0.d0
      !
      IF( j1 == j2 .AND. k1 == k2 ) Tnz(1) = .TRUE.
      IF( i1 == i2 .AND. k1 == k2 ) Tnz(2) = .TRUE.
      IF( i1 == i2 .AND. j1 == j2 ) Tnz(3) = .TRUE.
      !WRITE(*,*) ip1,ip2, any(Tnz)
      IF( ANY(Tnz) ) THEN
        NNZ = NNZ + 1
        IF( ip1 /= ip2 ) NNZ = NNZ + 1
      ENDIF
        !
    ENDDO
  ENDDO

  !WRITE(*,*) 'NNZ, Npoints = ', NNZ, Npoints
  !WRITE(*,*) 'Percentage   = ', ( NNZ ) / dble(Npoints**2) * 100d0
END FUNCTION
