!WRF:MODEL_LAYER:DYNAMICS
!
MODULE module_planetary_damping_em

  USE module_wrf_error

CONTAINS

!==============================================================================

  SUBROUTINE damp_field_at_top( case, field, field_tend,      &
                                mu, P2SI, t0, pi,             &
                                to_zonal_mean,                &
                                ids, ide, jds, jde, kds, kde, &
                                ims, ime, jms, jme, kms, kme, &
                                its, ite, jts, jte, kts, kte )

    !----------------------------------------------------------------
    ! Simulate gravity wave breaking by damping momentum and energy
    ! near the model top
    !----------------------------------------------------------------

    IMPLICIT NONE

    INTEGER,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte

    CHARACTER(LEN=1), INTENT(IN   ) :: case
    LOGICAL,          INTENT(IN   ) :: to_zonal_mean
    REAL,             INTENT(IN   ) :: P2SI, t0
    REAL, DIMENSION( ims:ime , jms:jme )           , INTENT(IN   ) ::  mu
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) ::  field, &
                                                                       pi
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) ::  field_tend

    ! Define some choices for what to do with model top damping
    ! to_zonal_mean = .TRUE.  if damping forces values to an instantaneous
    !                         zonal-mean value
    !               = .FALSE. if damping is forced to a pre-defined value
    !                         t_naught for T (which must be defined below)
    !                         0. for u and v
    !
    ! taudmp will be relaxation time in *MINUTES* ... we will take care of
    ! the conversion to 1/time (in units of 1/s) when making rtau(:)

    INTEGER, PARAMETER :: nlay_dmp = 3
    REAL, DIMENSION(nlay_dmp), PARAMETER :: taudmp = &     
          (/ 2880., 8640., 25920. /) ! 2 days, 6 days, 18 days
    REAL, PARAMETER :: t_naught = 140.

    REAL, DIMENSION(ims:ime, kms:kme, jms:jme) :: damped
    REAL, DIMENSION(nlay_dmp) :: rtau
    REAL :: rnx, t_damp, sum_damped

    INTEGER :: i, j, k, i_end, j_end
    INTEGER :: k_start, k_end

    IF (kte < kde-nlay_dmp) return

    k_start = kde-nlay_dmp
    k_end   = kde-1
    IF (kts > kde-nlay_dmp) k_start = kts
    IF (kte < kde-1       ) k_end   = kte

    rtau(:) = 1./(taudmp(:)*60.*P2SI)

    IF (to_zonal_mean .AND. ((its /= ids) .OR. (ite /= ide))) THEN
       WRITE ( wrf_err_message , * ) 'module_planetary_damping: damp_field_at_top: (its /= ids) or (ite /= ide)',its,ids,ite,ide
       CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
    END IF

    i_end = MIN(ite,ide-1)
    IF (case == 'U') i_end = ite
    j_end = MIN(jte,jde-1)
    IF (case == 'V') j_end = jte

    DO j=jts,j_end
    DO k=k_start,k_end
    DO i=ids,i_end
       damped(i,k,j)=0.
    END DO
    END DO
    END DO

    IF (case == 'U') THEN
       rnx = 1./REAL(ide-ids+1)
    ELSE IF ((case == 'V') .OR. (case == 'T')) THEN
       rnx = 1./REAL(ide-ids)
    ELSE
       WRITE ( wrf_err_message , * ) 'module_planetary_damping: ', &
            'damp_field_at_top: unknown field being damped: ',TRIM(case)
       CALL wrf_error_fatal ( TRIM( wrf_err_message ) )
    END IF
    IF (to_zonal_mean) THEN
       DO j=jts,j_end
       DO k=k_start,k_end
          sum_damped = 0.
          DO i=ids,i_end
             IF (case == 'T') THEN
                ! Damp to an average temperature
                sum_damped = sum_damped+(field(i,k,j)+t0)*pi(i,k,j)*rnx
             ELSE
                sum_damped = sum_damped+field(i,k,j)*rnx
             END IF
          END DO
          DO i=ids,i_end
             IF (case == 'T') THEN
                ! The field we are damping in case 'T' is actually the
                ! variable "t_2", the perturbation potential temperature,
                ! so convert the mean value back into a perturbation
                ! potential temperature
                damped(i,k,j) = (sum_damped/pi(i,k,j))-t0
             ELSE
                damped(i,k,j) = sum_damped
             END IF
          END DO
       END DO
       END DO
    ELSE
       ! We are modifying a potential temperature (t_2, ultimately)
       ! When we relax to an average, we are safe to average the potential
       ! temperatures
       ! When we relax to a constant temperature, "t_naught", we will need
       ! to convert that constant to a potential temperature
       IF (case == 'T') THEN
          DO j=jts,     j_end
          DO k=k_start, k_end
          DO i=ids,     i_end
             ! "damped" will be our damped *potential* temperature
             ! Need to convert t_naught to a perturbation potential temperature
             damped(i,k,j) = (t_naught / pi(i,k,j)) - t0
          END DO
          END DO
          END DO
       END IF
    END IF


    ! tendencies only need to be coupled with mu, but the correct mu variable
    ! for the different stagger points, which will be handled by passing in
    ! the correct array
    DO j = jts,     j_end
    DO k = k_start, k_end
    DO i = its,     i_end
       field_tend(i,k,j) =  field_tend(i,k,j) &
                           -(field(i,k,j)-damped(i,k,j))*rtau(kde-k)*mu(i,j)
    END DO
    END DO
    END DO

  END SUBROUTINE damp_field_at_top

!------------------------------------------------------------------------------

END MODULE module_planetary_damping_em
