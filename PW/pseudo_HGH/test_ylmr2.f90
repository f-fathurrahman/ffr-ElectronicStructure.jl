PROGRAM test_ylmr2

  IMPLICIT NONE
  INTEGER :: lmax2
  INTEGER :: Ng
  REAL(8), ALLOCATABLE :: G(:,:), GG(:)
  REAL(8), ALLOCATABLE :: ylm(:,:)
  INTEGER :: ig
  INTEGER :: lmaxkb, lm

  Ng = 5
  lmaxkb = 3
  lmax2 = (lmaxkb+1)**2

  WRITE(*,*) 'lmax2 = ', lmax2

  ALLOCATE( G(3,Ng) )
  ALLOCATE( GG(Ng) )
  ALLOCATE( ylm(Ng,lmax2) )

  DO ig = 1, Ng
    G(:,ig) = (/ -2.1d0, 2.d0, 1.d0 /)
    GG(ig) = G(1,ig)**2 + G(2,ig)**2 + G(3,ig)**2
  ENDDO

  CALL ylmr2( lmax2, Ng, G, GG, ylm )

  DO lm = 1,lmax2
    WRITE(*,*)
    WRITE(*,'(A,I5)') 'lm = ', lm
    DO ig = 1,Ng
      WRITE(*,'(I5,F18.10)') ig, ylm(ig,lm)
    ENDDO
  ENDDO

  DEALLOCATE( G )
  DEALLOCATE( GG )
  DEALLOCATE( ylm )

END PROGRAM
