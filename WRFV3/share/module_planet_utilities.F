MODULE module_planet_utilities

  USE module_model_constants

  PRIVATE
  PUBLIC :: get_ls

  REAL(KIND(0d0)), PARAMETER :: pi_d = 3.141592653589793238462643383279502884197d0

CONTAINS

!----------------------------------------------------------------
  REAL FUNCTION get_ls(julian) RESULT (ls)
!----------------------------------------------------------------
!
! Calculate solar longitude as a function of Julian day.
! The definition used for this subroutine (and similarly for
! the calculation in phys/module_radation_driver.F:radconst is
! midnight at the beginning of the first day of the year is
! equal to a "Julian day" of 0.000...
!
! Input: Julian day (sols, fractional)
!
! Output: Solar longitude (degrees)
!
!----------------------------------------------------------------
    IMPLICIT NONE
!----------------------------------------------------------------

!   Input/Ouptut variables
    REAL, INTENT(IN) :: julian     ! Day of the year, sols

!   Parameter variables
    REAL(KIND(0d0)), PARAMETER :: small_value = 1.D-6

!   Local Variables
    REAL(KIND(0d0)) :: deleqn, date_dbl
    REAL(KIND(0d0)) :: er, qq, e, cd0, ep, em
    REAL(KIND(0d0)) :: eq, w, als

    deleqn = equinox_fraction * REAL(planet_year) 

    !-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:
        
    !  DATE = DAYS SINCE LAST PERIHELION PASSAGE
    date_dbl = REAL(julian,KIND(0.d0)) - zero_date
    DO WHILE (date_dbl < 0.)
       date_dbl=date_dbl+REAL(planet_year)
    END DO
    DO WHILE (date_dbl > REAL(planet_year))
       date_dbl=date_dbl-REAL(planet_year)
    END DO

    er = SQRT( (1.d0+eccentricity)/(1.d0-eccentricity) )

    !  qq is the mean anomaly
    qq = 2.d0 * (pi_d * deleqn / REAL(planet_year))

    !  determine true anomaly at equinox:  eq
    !  Iteration for eq
    e = 1.d0
    cd0 = 1.d0
    DO WHILE (cd0 > small_value)
       ep = e - (e-eccentricity*SIN(e)-qq)/(1.d0-eccentricity*COS(e))
       cd0 = ABS(e-ep)
       e = ep
    END DO
    eq = 2.d0 * ATAN( er * TAN(0.5d0*e) )

    !  determine true anomaly at current date:  w
    !  Iteration for w
    em = 2.d0 * pi_d * date_dbl / REAL(planet_year)
    e = 1.d0
    cd0 = 1.d0
    DO WHILE (cd0 > small_value)
       ep = e-(e-eccentricity*SIN(e)-em)/(1.d0-eccentricity*COS(e))
       cd0 = ABS(e-ep)
       e = ep
    END DO
    w = 2.d0 * ATAN( er * TAN(0.5d0*e) )

    !  Radius vector ( astronomical units:  AU )
    als= (w - eq)*180.d0/pi_d      ! Solar Longitude
    IF (als < 0.d0) als=als+360.d0
    ls = REAL(als)
  
  END FUNCTION get_ls

END MODULE module_planet_utilities
