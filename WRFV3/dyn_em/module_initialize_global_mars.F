!IDEAL:MODEL_LAYER:INITIALIZATION

!  This MODULE holds the routines which are used to perform various initializations
!  for the individual domains.  

!-----------------------------------------------------------------------

MODULE module_initialize_ideal

   USE module_domain
   USE module_io_domain
   USE module_model_constants
   USE module_bc
   USE module_wrf_error
   USE module_planet_utilities
#ifdef DM_PARALLEL
   USE module_dm
#endif


CONTAINS


!-------------------------------------------------------------------
! this is a wrapper for the solver-specific init_domain routines.
! Also dereferences the grid variables and passes them down as arguments.
! This is crucial, since the lower level routines may do message passing
! and this will get fouled up on machines that insist on passing down
! copies of assumed-shape arrays (by passing down as arguments, the 
! data are treated as assumed-size -- ie. f77 -- arrays and the copying
! business is avoided).  Fie on the F90 designers.  Fie and a pox.

   SUBROUTINE init_domain ( grid )

   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid 
   !  Local data.
   INTEGER :: idum1, idum2

   CALL set_scalar_indices_from_config ( head_grid%id , idum1, idum2 )

     CALL init_domain_rk( grid &
!
#include <actual_new_args.inc>
!
                        )

   END SUBROUTINE init_domain

!-------------------------------------------------------------------

   SUBROUTINE init_domain_rk ( grid &
!
# include <dummy_new_args.inc>
!
)
   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid

# include <dummy_decl.inc>

   TYPE (grid_config_rec_type)              :: config_flags

   !  Local data
   INTEGER                             ::                       &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte, &
                                  i, j, k, l

   INTEGER ::              ips , ipe , jps , jpe , kps , kpe

   INTEGER ::              imsx, imex, jmsx, jmex, kmsx, kmex,    &
                           ipsx, ipex, jpsx, jpex, kpsx, kpex,    &
                           imsy, imey, jmsy, jmey, kmsy, kmey,    &
                           ipsy, ipey, jpsy, jpey, kpsy, kpey

   ! Local data

   ! Some character strings to report setup conditions and progress
   CHARACTER(LEN=80) :: description
   CHARACTER(LEN=80) :: string1

   INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
   INTEGER :: nxx, nyy, error, icm, jcm, idealloc
   INTEGER :: ig, jg
   REAL :: icent, jcent
   REAL :: dphi, dlam
   REAL :: vlat, ulong
   REAL :: p_surf, p_level, pd_surf, qvf1, qvf2, qvf, tperturb
   REAL :: thtmp, ptmp, temp(3), cof1, cof2

   INTEGER :: time_step, time_step_fract_num, time_step_fract_den
   REAL :: dt


   CALL get_ijk_from_grid (  grid ,                   &
                             ids, ide, jds, jde, kds, kde,    &
                             ims, ime, jms, jme, kms, kme,    &
                             ips, ipe, jps, jpe, kps, kpe,    &
                             imsx, imex, jmsx, jmex, kmsx, kmex,    &
                             ipsx, ipex, jpsx, jpex, kpsx, kpex,    &
                             imsy, imey, jmsy, jmey, kmsy, kmey,    &
                             ipsy, ipey, jpsy, jpey, kpsy, kpey )

   SELECT CASE ( model_data_order )
         CASE ( DATA_ORDER_ZXY )
   kds = grid%sd31 ; kde = grid%ed31 ;
   ids = grid%sd32 ; ide = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   kms = grid%sm31 ; kme = grid%em31 ;
   ims = grid%sm32 ; ime = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   kts = grid%sp31 ; kte = grid%ep31 ;   ! note that tile is entire patch
   its = grid%sp32 ; ite = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XYZ )
   ids = grid%sd31 ; ide = grid%ed31 ;
   jds = grid%sd32 ; jde = grid%ed32 ;
   kds = grid%sd33 ; kde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   jms = grid%sm32 ; jme = grid%em32 ;
   kms = grid%sm33 ; kme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   jts = grid%sp32 ; jte = grid%ep32 ;   ! note that tile is entire patch
   kts = grid%sp33 ; kte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XZY )
   ids = grid%sd31 ; ide = grid%ed31 ;
   kds = grid%sd32 ; kde = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   kms = grid%sm32 ; kme = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   kts = grid%sp32 ; kte = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch

   END SELECT

   CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )

! here we check to see if the boundary conditions are set properly

   CALL boundary_condition_check( config_flags, bdyzone, error, grid%id )

   grid%itimestep=0
   grid%step_number = 0

   CALL nl_get_time_step ( grid%id, time_step )
   CALL nl_get_time_step_fract_num( grid%id, time_step_fract_num )
   CALL nl_get_time_step_fract_den( grid%id, time_step_fract_den )
   dt = REAL(time_step) + REAL(time_step_fract_num) / REAL(time_step_fract_den)
   dt = dt * P2SI
   CALL nl_set_dt( grid%id, dt )
   grid%dt = dt

#ifdef DM_PARALLEL
   CALL wrf_dm_bcast_bytes( icm , IWORDSIZE )
   CALL wrf_dm_bcast_bytes( jcm , IWORDSIZE )
#endif

   ! -------------------------------------------------------------------
   ! Set horizontal grid basic definitions

   nyy = jde-jds
   nxx = ide-ids
   icent = REAL(nxx+1)/2.
   jcent = REAL(nyy+1)/2.
   dphi = config_flags%lat_span/REAL(nyy)
   dlam = config_flags%lon_span/REAL(nxx)
   DO j = jts, jte
      DO i = its, ite
         ! ig is the I index in the global (domain) span of the array.
         ! jg is the J index in the global (domain) span of the array.
         ig = i - ids + 1
         jg = j - jds + 1
         grid%xlat(i,j)  = (REAL(jg)-jcent)*dphi+config_flags%lat_ref
         grid%xlong(i,j) = (REAL(ig)-icent)*dlam+config_flags%lon_ref
         vlat       = grid%xlat(i,j)  - 0.5*dphi
         ulong      = grid%xlong(i,j) - 0.5*dlam
         DO
            IF (grid%xlong(i,j) >= -180.) EXIT
            grid%xlong(i,j) = grid%xlong(i,j) + 360.
         END DO
         DO
            IF (grid%xlong(i,j) <=  180.) EXIT
            grid%xlong(i,j) = grid%xlong(i,j) - 360.
         END DO
         DO
            IF (ulong >= -180.) EXIT
            ulong = ulong + 360.
         END DO
         DO
            IF (ulong <=  180.) EXIT
            ulong = ulong - 360.
         END DO
         grid%xlat_u(i,j)  = grid%xlat(i,j)
         grid%xlong_u(i,j) = ulong
         grid%xlat_v(i,j)  = vlat
         grid%xlong_v(i,j) = grid%xlong(i,j)

         grid%msftx(i,j) = 1./COS(grid%xlat(i,j)*degrad)
         grid%msfty(i,j) = 1.
         grid%msfux(i,j) = 1./COS(grid%xlat(i,j)*degrad)
         grid%msfuy(i,j) = 1.
         grid%msfvx_inv(i,j) = COS(vlat*degrad)
         grid%msfvx(i,j) = 1./grid%msfvx_inv(i,j)
         grid%msfvy(i,j) = 1.
         IF (config_flags%lat_span == 180.) THEN
            IF ( (jts == jds) .AND. (j == jts) )THEN
               grid%msfvx(i,j)     = 0.
               grid%msfvx_inv(i,j) = 0.
            END IF
            IF ( (jte == jde) .AND. (j == jte) )THEN
               grid%msfvx(i,j)     = 0.
               grid%msfvx_inv(i,j) = 0.
            END IF
         END IF

         grid%e(i,j)     = 2*EOMEG*COS(grid%xlat(i,j)*degrad)
         grid%f(i,j)     = 2*EOMEG*SIN(grid%xlat(i,j)*degrad)

         grid%sina(i,j) = 0.
         grid%cosa(i,j) = 1.

         grid%clat(i,j)  = grid%xlat(i,j)
         grid%clong(i,j) = grid%xlong(i,j)

      END DO
   END DO

   ! -------------------------------------------------------------------
   ! Generate L_s for output to file
   grid%julian = config_flags%start_day - 1       + &
                (config_flags%start_hour/24.)     + &
                (config_flags%start_minute/1440.) + &
                (config_flags%start_second/86400.)
   grid%l_s = GET_LS(grid%julian)
   description = 'L_s:'
   WRITE(string1,'(F16.8)') grid%l_s
   WRITE(wrf_err_message,*) TRIM(description),'                 ',TRIM(string1)
   CALL wrf_message ( TRIM( wrf_err_message ) )

   ! -------------------------------------------------------------------
   ! Namelist definitions

   grid%rdx = 1./grid%dx
   grid%rdy = 1./grid%dy

   description = 'DX:'
   WRITE(string1,'(F16.8)') grid%dx
   WRITE(wrf_err_message,*) TRIM(description),'                  ', &
                            TRIM(string1)
   CALL wrf_message ( TRIM( wrf_err_message ) )
   description = 'DY:'
   WRITE(string1,'(F16.8)') grid%dy
   WRITE(wrf_err_message,*) TRIM(description),'                  ', &
                            TRIM(string1)
   CALL wrf_message ( TRIM( wrf_err_message ) )
   WRITE( wrf_err_message , * ) ''
   CALL wrf_message ( TRIM( wrf_err_message ) )


   grid%cen_lat      = config_flags%lat_ref
   grid%cen_lon      = config_flags%lon_ref
   grid%truelat1     = 0.
   grid%truelat2     = 0.
   grid%moad_cen_lat = config_flags%lat_ref
   grid%stand_lon    = config_flags%lon_ref
   grid%map_proj     = 0
   grid%iswater      = 0
   CALL nl_set_cen_lat ( grid%id , grid%cen_lat )
   CALL nl_set_cen_lon ( grid%id , grid%cen_lon )
   CALL nl_set_truelat1 ( grid%id , grid%truelat1 )
   CALL nl_set_truelat2 ( grid%id , grid%truelat2 )
   CALL nl_set_moad_cen_lat ( grid%id , grid%moad_cen_lat )
   CALL nl_set_stand_lon ( grid%id , grid%stand_lon )
   CALL nl_set_map_proj ( grid%id , grid%map_proj )
   ! -------------------------------------------------------------------


   ! -------------------------------------------------------------------
   ! Do the 2D surface arrays
   !
   DO j=jts,jte
   DO i=its,ite
      grid%xland(i,j)    = 1.  ! 1 for land, 2 for water
      grid%landmask(i,j) = 1.  ! 1 for land, 0 for water
      grid%mavail(i,j)   = 0.
      grid%angslope(i,j) = 0.
      grid%azmslope(i,j) = 0.
      grid%ht(i,j)       = 0.
      grid%albedo(i,j)   = 0.2
      grid%thc(i,j)      = 216.
      grid%znt(i,j)      = 0.023
      grid%emiss(i,j)    = 1.
   END DO
   END DO

   ! -------------------------------------------------------------------
   ! Set up the vertical grid definitions
   DO k = kds, kde
      grid%znw(k) = 1. - REAL(k-kds)/REAL(kde-kds)
   END DO
   DO k=1, kde-1
      grid%dnw(k) = grid%znw(k+1) - grid%znw(k)
      grid%rdnw(k) = 1./grid%dnw(k)
      grid%znu(k) = 0.5*(grid%znw(k+1)+grid%znw(k))
   END DO
   DO k=2, kde-1
      grid%dn(k) = 0.5*(grid%dnw(k)+grid%dnw(k-1))
      grid%rdn(k) = 1./grid%dn(k)
      grid%fnp(k) = .5* grid%dnw(k  )/grid%dn(k)
      grid%fnm(k) = .5* grid%dnw(k-1)/grid%dn(k)
   END DO

   cof1 = (2.*grid%dn(2)+grid%dn(3))/(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(2) 
   cof2 =     grid%dn(2)        /(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(3) 
   grid%cf1  = grid%fnp(2) + cof1
   grid%cf2  = grid%fnm(2) - cof1 - cof2
   grid%cf3  = cof2       

   grid%cfn  = (.5*grid%dnw(kde-1)+grid%dn(kde-1))/grid%dn(kde-1)
   grid%cfn1 = -.5*grid%dnw(kde-1)/grid%dn(kde-1)


   ! -------------------------------------------------------------------
   ! Set values needed when spawning the domain and for use in
   ! calculations below
   grid%t00  = config_flags%base_temp
   grid%p00  = config_flags%base_pres
   grid%tlp  = config_flags%base_lapse
   grid%tiso = config_flags%iso_temp

   ! -------------------------------------------------------------------
   ! find ptop for the desired ztop (ztop is input from the namelist)
   grid%p_top = Z2P(config_flags%ztop, config_flags%base_pres,       &
                    config_flags%base_temp, config_flags%base_lapse, &
                    config_flags%iso_temp, g, R_d)

   i_start = its
   i_end   = ite
   j_start = jts
   j_end   = jte

   ! Values of geopotential (base, perturbation, and at p0) at the surface
   DO j=j_start,j_end
   DO i=i_start,i_end
      grid%phb(i,1,j) = grid%ht(i,j)*g
      grid%php(i,1,j) = 0.
      grid%ph0(i,1,j) = grid%ht(i,j)*g
   ENDDO
   ENDDO


   ! Calculate base state total mass (mub), pressure (pb), inverse
   ! density (alb = 1/rho), and geopotential (phb).
   DO J = j_start, j_end
   DO I = i_start, i_end

      p_surf = Z2P(grid%phb(i,1,j)/g, config_flags%base_pres,       &
                   config_flags%base_temp, config_flags%base_lapse, &
                   config_flags%iso_temp, g, R_d)
      grid%mub(i,j) = p_surf-grid%p_top

      ! given p (coordinate), calculate theta and compute 1/rho from equation
      ! of state

      DO K = kts, kte-1
         p_level = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
         grid%pb(i,k,j) = p_level
         grid%t_init(i,k,j) = P2T(p_level, config_flags%base_pres, &
                                   config_flags%base_temp,         &
                                   config_flags%base_lapse,        &
                                   config_flags%iso_temp, g, R_d)
         grid%t_init(i,k,j) = grid%t_init(i,k,j)*(p0/p_level)**rcp
         ! Remember the all important t0!
         grid%t_init(i,k,j) = grid%t_init(i,k,j) - t0

         ! Assuming p0 = p1000mb...
         grid%alb(i,k,j)=(r_d/p1000mb)*(grid%t_init(i,k,j)+t0)*(grid%pb(i,k,j)/p1000mb)**cvpm
      END DO

      ! calculate hydrostatic balance (alternatively we could interpolate
      ! the geopotential from the sounding, but this assures that the base
      ! state is in exact hydrostatic balance with respect to the model eqns.)

      DO k = kts+1, kte
         grid%phb(i,k,j) = grid%phb(i,k-1,j) - grid%dnw(k-1)*grid%mub(i,j)*grid%alb(i,k-1,j)
      ENDDO

   ENDDO
   ENDDO

   ! Initialize moisture (water vapor)
   DO j = j_start, j_end
      DO k = kts, kte
         DO i = i_start, i_end
            moist(i,k,j,P_QV) = 0.
         END DO
      END DO
   END DO


   ! Now calculate the full (hydrostatically-balanced) state for each column
   ! We will include moisture
   DO J = j_start, j_end
   DO I = i_start, i_end

      ! At this point p_top is already set. find the DRY mass in the column
      ! Use the z->p function to find the surface pressure from the surface
      ! topography value (z=ht(i,j)=phb(i,1,j)/g)
      pd_surf = Z2P(grid%phb(i,1,j)/g, config_flags%base_pres,       &
                    config_flags%base_temp, config_flags%base_lapse, &
                    config_flags%iso_temp, g, R_d)
      grid%psfc(i,j) = pd_surf

      ! compute the perturbation mass (mu/mu_1/mu_2) and the full mass
      grid%mu_1(i,j) = pd_surf-grid%p_top - grid%mub(i,j)
      ! The following was used to put an initial 'ping' into the system,
      ! much like dropping a stone into a pond, to check e.g. the polar BCs
      ! IF((i == ide/2) .and. (j == jde/2)) mu_1(i,j)=mub(i,j)*.05
      grid%mu_2(i,j) = grid%mu_1(i,j)
      grid%mu0(i,j)  = grid%mu_1(i,j) + grid%mub(i,j)

      ! given the dry pressure and coordinate system, calculate the
      ! perturbation potential temperature (t/t_1/t_2)

      DO k = kds, kde-1
         p_level = grid%znu(k)*(pd_surf - grid%p_top) + grid%p_top
         ! Use the p->T function to find the temperature given a pressure
         grid%t_1(i,k,j) = P2T(p_level, config_flags%base_pres,                 &
                               config_flags%base_temp, config_flags%base_lapse, &
                               config_flags%iso_temp, g, R_d)
         grid%t_1(i,k,j) = grid%t_1(i,k,j)*(p0/p_level)**rcp
         grid%t_1(i,k,j) = grid%t_1(i,k,j)-t0
         grid%t_2(i,k,j) = grid%t_1(i,k,j)
      END DO


      ! integrate the hydrostatic equation (from the RHS of the bigstep
      ! vertical momentum equation) down from the top to get p.
      ! first from the top of the model to the top pressure

      k = kte-1  ! top level

      qvf1 = 0.5*(grid%moist(i,k,j,P_QV)+grid%moist(i,k,j,P_QV))
      qvf2 = 1./(1.+qvf1)
      qvf1 = qvf1*qvf2

      grid%p(i,k,j) = - 0.5*(grid%mu_1(i,j)+qvf1*grid%mub(i,j))/grid%rdnw(k)/qvf2
      qvf = 1. + rvovrd*grid%moist(i,k,j,P_QV)
      grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                  (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
      grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)

      !  down the column

      do k=kte-2,kts,-1
         qvf1 = 0.5*(grid%moist(i,k,j,P_QV)+grid%moist(i,k+1,j,P_QV))
         qvf2 = 1./(1.+qvf1)
         qvf1 = qvf1*qvf2
         grid%p(i,k,j) = grid%p(i,k+1,j) - (grid%mu_1(i,j) + qvf1*grid%mub(i,j))/qvf2/grid%rdn(k+1)
         qvf = 1. + rvovrd*grid%moist(i,k,j,P_QV)
         grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                     (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
         grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
      enddo

      ! this is the hydrostatic equation used in the model after the
      ! small timesteps.  In the model, al (inverse density)
      ! is computed from the geopotential.

      grid%ph_1(i,1,j) = 0.
      DO k  = kts+1,kte
         grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (1./grid%rdnw(k-1))*(  &
                      (grid%mub(i,j)+grid%mu_1(i,j))*grid%al(i,k-1,j)+  &
                       grid%mu_1(i,j)*grid%alb(i,k-1,j)  )

         grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
         grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

   END DO
   END DO

   ! Initialize winds
   DO j=j_start,j_end
   DO k=kds,kde
   DO i=i_start,i_end
      grid%u_1(i,k,j) = 0.
      grid%u_2(i,k,j) = 0.
      grid%v_1(i,k,j) = 0.
      grid%v_2(i,k,j) = 0.
      grid%ww(i,k,j)  = 0.
      grid%w_1(i,k,j) = 0.
      grid%w_2(i,k,j) = 0.
      grid%h_diabatic(i,k,j) = 0.
   END DO
   END DO
   END DO

  
   ! -------------------------------------------------------------------
   ! The "base" variables.
   DO k=kts,kte
      grid%t_base(k)  = 0.
      grid%qv_base(k) = 0.
      grid%u_base(k)  = 0.
      grid%v_base(k)  = 0.
   END DO

   ! -------------------------------------------------------------------
   ! Do the subsurface arrays

   ! Only one layer: infinite slab at constant temperature below the surface
   ! Surface temperature is an infinitely thin "skin" on top of an
   ! half-infinite slab.
   ! And the temperature of both the skin and the slab are determined from
   ! the initial nearest-surface-air-layer temperature.
   DO J = jts, MIN(jte, jde-1)
      DO I = its, MIN(ite, ide-1)
         thtmp   = grid%t_2(i,1,j)+t0
         ptmp    = grid%p(i,1,j)+grid%pb(i,1,j)
         temp(1) = thtmp * (ptmp/p1000mb)**rcp
         thtmp   = grid%t_2(i,2,j)+t0
         ptmp    = grid%p(i,2,j)+grid%pb(i,2,j)
         temp(2) = thtmp * (ptmp/p1000mb)**rcp
         thtmp   = grid%t_2(i,3,j)+t0
         ptmp    = grid%p(i,3,j)+grid%pb(i,3,j)
         temp(3) = thtmp * (ptmp/p1000mb)**rcp
         grid%tsk(I,J)=cf1*temp(1)+cf2*temp(2)+cf3*temp(3)
         grid%tmn(I,J)=grid%tsk(I,J)-0.5
      END DO
   END DO

   WRITE( wrf_err_message , * ) ''
   CALL wrf_message ( TRIM( wrf_err_message ) )

   grid%save_topo_from_real=1

   RETURN

 END SUBROUTINE init_domain_rk

!---------------------------------------------------------------------

 SUBROUTINE init_module_initialize
 END SUBROUTINE init_module_initialize

!---------------------------------------------------------------------

 FUNCTION Z2P(z, base_pres, base_temp, base_lapse, iso_temp, g, R_d) RESULT (p)
   IMPLICIT NONE
   REAL, INTENT(IN) :: z
   REAL, INTENT(IN) :: base_pres, base_temp, base_lapse, iso_temp, g, R_d
   REAL             :: p
   REAL             :: z_iso, p_iso
   z_iso = R_d*(base_temp*base_temp - iso_temp*iso_temp)/(2.*g*base_lapse)
   IF (z < z_iso) THEN
      p = base_pres*EXP(-(base_temp/base_lapse) * &
                         (1.-SQRT(1.-2.*g*base_lapse*z / &
                                  (R_d*base_temp*base_temp))))
   ELSE
      p_iso = base_pres*EXP(-(base_temp/base_lapse) * &
                             (1.-SQRT(1.-2.*g*base_lapse*z_iso / &
                                      (R_d*base_temp*base_temp))))
      p     = p_iso * EXP(-(g*(z-z_iso))/(R_d*iso_temp))
   END IF
 END FUNCTION Z2P

 FUNCTION P2T(p, base_pres, base_temp, base_lapse, iso_temp, g, R_d) RESULT (t)
   REAL, INTENT(IN) :: p
   REAL, INTENT(IN) :: base_pres, base_temp, base_lapse, iso_temp, g, R_d
   REAL             :: t
   t = base_temp + base_lapse*LOG(p/base_pres)
 END FUNCTION P2T

END MODULE module_initialize_ideal
