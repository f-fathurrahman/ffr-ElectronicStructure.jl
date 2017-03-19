PROGRAM test_ylmr2

  IMPLICIT NONE 
  INTEGER :: lmax2
  INTEGER :: Ng
  REAL(8), ALLOCATABLE :: G(:,:), GG(:)
  REAL(8), ALLOCATABLE :: ylm(:,:)
  INTEGER :: ig
  INTEGER :: lmaxkb

  Ng = 5
  lmaxkb = 2
  lmax2 = (lmaxkb+1)**2

  WRITE(*,*) 'lmax2 = ', lmax2

  ALLOCATE( G(3,Ng) )
  ALLOCATE( GG(Ng) )
  ALLOCATE( ylm(Ng,lmax2) )

  DO ig = 1, Ng
    G(:,ig) = (/ 1.3d0*ig, 1.1d0*ig, 2.1d0*ig /)
    GG(ig) = G(1,ig)**2 + G(2,ig)**2 + G(3,ig)**2
  ENDDO 
  
  CALL ylmr2( lmax2, Ng, G, GG, ylm )

  DEALLOCATE( G )
  DEALLOCATE( GG )
  DEALLOCATE( ylm )

END PROGRAM 

