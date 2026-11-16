! Dummy module for cloud fraction scheme to ensure that MONC and KiD build

MODULE cloud_frac_scheme

IMPLICIT NONE

CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='CLOUD_FRAC_SCHEME'

CONTAINS

SUBROUTINE cloud_frac_casim_mphys(k, pressure, T_in, T_l, rhcrit_lev, &
                                  qs, qv, cloud_mass, qfields, cloud_mass_new)

USE variable_precision, ONLY: wp

USE vernier_mod, ONLY: vernier_start, vernier_stop, vik
USE vernier_timing_mod, ONLY: LPROF_CASIM
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER,  INTENT(IN)  :: k
REAL(wp), INTENT(IN)  :: pressure, T_in, T_l, rhcrit_lev, qs, qv, cloud_mass, qfields
REAL(wp), INTENT(OUT) :: cloud_mass_new

CHARACTER(LEN=*),   PARAMETER :: RoutineName='CLOUD_FRAC_CASIM_MPHYS'
INTEGER(KIND=vik) :: vnr_handle_1

!--------------------------------------------------------------------------
! End of header, no more declarations beyond here
!--------------------------------------------------------------------------
IF (LPROF_CASIM) CALL vernier_start(vnr_handle_1, ModuleName//':'//RoutineName)

cloud_mass_new = cloud_mass

IF (LPROF_CASIM) CALL vernier_stop(vnr_handle_1)

END SUBROUTINE cloud_frac_casim_mphys

END MODULE cloud_frac_scheme
