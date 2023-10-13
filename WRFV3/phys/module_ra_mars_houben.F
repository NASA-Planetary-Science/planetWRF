!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_ra_mars_houben

! This module provides a very simple, newtonian cooling scheme for a Mars GCM
! It reads pre-calculated radiative convective model temperatures and yields
! a simulation with nearly identical forcing to that described by Haberle et
! al., 1997 (JGR-planets), and Houben et al., 1997 (JGR-planets)

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: houbeninit, houben

  INTEGER, PARAMETER ::  hh_nt = 671
  INTEGER, PARAMETER ::  hh_nz = 16
  INTEGER, PARAMETER ::  hh_nl = 32
  REAL, SAVE, DIMENSION(hh_nz,hh_nl,hh_nt) ::  hh_teqm
  REAL, SAVE, DIMENSION(hh_nz)             ::  hh_z
  REAL, SAVE, DIMENSION(hh_nl)             ::  hh_lats
  REAL, SAVE, DIMENSION(hh_nt)             ::  hh_ls

CONTAINS

!--------------------------------------------------------------------
  SUBROUTINE houbeninit(RTHRATEN,RTHRATENSW,RTHRATENLW,restart,  &
                        ids, ide, jds, jde, kds, kde,            &
                        ims, ime, jms, jme, kms, kme,            &
                        its, ite, jts, jte, kts, kte             )

    IMPLICIT NONE
    LOGICAL , INTENT(IN)           :: restart
    INTEGER , INTENT(IN)           :: ids, ide, jds, jde, kds, kde,  &
                                      ims, ime, jms, jme, kms, kme,  &
                                      its, ite, jts, jte, kts, kte

    REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) ::   &
                                                             RTHRATEN,   &
                                                             RTHRATENLW, &
                                                             RTHRATENSW
    INTEGER, PARAMETER :: unit = 2001
    INTEGER :: i, j, k, itf, jtf, ktf, itime, ilat

    jtf=min0(jte,jde-1)
    ktf=min0(kte,kde-1)
    itf=min0(ite,ide-1)

    IF(.NOT.restart)THEN
       DO j=jts,jtf
       DO k=kts,ktf
       DO i=its,itf
          RTHRATEN(i,k,j)=0.
          RTHRATENSW(i,k,j)=0.
          RTHRATENLW(i,k,j)=0.
       END DO
       END DO
       END DO
    END IF

    OPEN(unit,file='./Data/houben.txt',status='old',action='read')
    DO itime=2,hh_nt-1
       READ(unit,'(f8.2)') hh_ls(itime)
       READ(unit,'(16f8.3)') hh_z
       READ(unit,'(32f8.2)') hh_lats
       DO ilat=1,hh_nl
          READ(unit,'(16f6.1)') hh_teqm(:,ilat,itime)
       END DO
    ENDDO
    CLOSE(unit)

    ! Add padding to east/west edges
    hh_teqm(:,:,1)     = hh_teqm(:,:,hh_nt-1)
    hh_teqm(:,:,hh_nt) = hh_teqm(:,:,2)

    hh_ls(1)    =hh_ls(hh_nt-1)-360.
    hh_ls(hh_nt)=hh_ls(2)+360.

    ! Convert to meters
    hh_z(:)=hh_z(:)*1000.

  END SUBROUTINE houbeninit

!------------------------------------------------------------------------------
  SUBROUTINE houben( RTHRATEN, t, z, z_at_w, pi3d, &
                     dt, xlat, degrad, l_s,        &
                     ids, ide, jds, jde, kds, kde, &
                     ims, ime, jms, jme, kms, kme, &
                     its, ite, jts, jte, kts, kte)

    IMPLICIT NONE

    INTEGER ,       INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte

    REAL, DIMENSION( ims:ime , kms:kme, jms:jme ) , INTENT(INOUT) ::  RTHRATEN
    REAL, DIMENSION( ims:ime , kms:kme, jms:jme ) , INTENT(IN   ) ::  z, t,   &
                                                                      z_at_w, &
                                                                      pi3d
    REAL, DIMENSION( ims:ime , jms:jme )          , INTENT(IN   ) ::  xlat
    REAL,                                           INTENT(IN   ) ::  dt, l_s, &
                                                                      degrad
    REAL, DIMENSION( kms:kme ) :: znu


    REAL, DIMENSION( hh_nz , hh_nl )        :: teqm1, teqm2
    REAL, PARAMETER :: damp_time = 2.*86400. ! 2 sols -> seconds
    REAL, PARAMETER :: vfac = 1./damp_time

    INTEGER :: i, j, k, it, il, iz
    REAL    :: t_naught, t01, t02, tfrac1, tfrac2

    ! Find time index ("it") such that hh_ls(it) <= l_s <= hh_ls(it+1)
    it = LOCATE(hh_ls,l_s)

    ! File is already padded with wrap arounds, so we are safe
    ! in limiting out of bounds back to edges
    IF (it < 1)       it=1
    IF (it > hh_nt-1) it=hh_nt-1
    tfrac1 = (l_s-hh_ls(it+1))/(hh_ls(it)-hh_ls(it+1))
    tfrac2 = (l_s-hh_ls(it  ))/(hh_ls(it+1)-hh_ls(it))

    teqm1 = hh_teqm(:,:,it)
    teqm2 = hh_teqm(:,:,it+1)

    DO j=jts,jte
    DO k=kts,kte
    DO i=its,ite
       t01 = BILINEAR(hh_z,hh_lats,teqm1,z(i,k,j),xlat(i,j))
       t02 = BILINEAR(hh_z,hh_lats,teqm2,z(i,k,j),xlat(i,j))
       t_naught = tfrac1*t01 + tfrac2*t02
       RTHRATEN(i,k,j)=RTHRATEN(i,k,j)-vfac*(t(i,k,j)-t_naught)/pi3D(i,k,j)
    END DO
    END DO
    END DO

  END SUBROUTINE houben

!--------------------------------------------------------------------
  REAL FUNCTION bilinear(xa,ya,za,x,y)
    ! Performs bilinear interpolation
    IMPLICIT NONE

    ! Input/output variables
    REAL,    INTENT(IN   )                 :: x,y
    REAL,    INTENT(IN   ), DIMENSION(:)   :: xa, ya
    REAL,    INTENT(IN   ), DIMENSION(:,:) :: za

    ! Local variables
    INTEGER :: m,n               ! Array sizes
    INTEGER :: j,k               ! Indices of locations near x,y
    REAL    :: t,u,z1,z2,z3,z4   ! Temporary values

    ! Find M, N, the size of the input arrays
    m = SIZE(xa)
    n = SIZE(ya)

    ! Find j and k such that xa(j) <= x <= xa(j+1) and
    !                        ya(k) <= y <= ya(k+1)
    j = LOCATE(xa,x)
    k = LOCATE(ya,y)

    ! Here we are saying if our value lies outside our box range,
    ! just use the value at the edge of the range
    IF (j < 1)   j=1
    IF (j > m-1) j=m-1
    IF (k < 1)   k=1
    IF (k > n-1) k=n-1

    ! Define the values of the four points around X,Y
    z1=za(j  ,k  )
    z2=za(j+1,k  )
    z3=za(j+1,k+1)
    z4=za(j  ,k+1)

    ! Set up variables for convenience
    t =(x-xa(j))/(xa(j+1)-xa(j))
    u =(y-ya(k))/(ya(k+1)-ya(k))
    IF (t < 0.) t = 0.
    IF (t > 1.) t = 1.
    IF (u < 0.) u = 0.
    IF (u > 1.) u = 1.

    ! Calculate the final value
    bilinear = (1.-t)*(1.-u)*z1 + &
                   t *(1.-u)*z2 + &
                   t *    u *z3 + &
               (1.-t)*    u *z4

  END FUNCTION bilinear

!--------------------------------------------------------------------
  INTEGER FUNCTION locate(y,x)
    IMPLICIT NONE
    ! Input/output variables
    REAL, INTENT(IN), DIMENSION(:) :: y
    REAL, INTENT(IN) :: x
    ! Local variables
    INTEGER :: n,jl,jm,ju
    LOGICAL :: ascnd

    n=SIZE(y)
    ascnd = (y(n) >= y(1))
    jl=0
    ju=n+1
    DO
       IF (ju-jl <= 1) EXIT
       jm=(ju+jl)/2
       IF (ascnd .EQV. (x >= y(jm))) THEN
          jl=jm
       ELSE
          ju=jm
       END IF
    END DO
    IF (x == y(1)) THEN
       locate=1
    ELSE IF (x == y(n)) THEN
       locate=n-1
    ELSE
       locate=jl
    END IF
  END FUNCTION locate

END MODULE module_ra_mars_houben
