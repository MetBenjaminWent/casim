module ventfac
  use variable_precision, only: wp
  use passive_fields, only: rho
  use special, only: pi, GammaFunc
  use mphys_parameters, only: hydro_params, vent_1, vent_2
  use mphys_constants, only: visair, rho0, Dv

!prf
  use mphys_switches, only: l_kfsm
!prf


  implicit none
  private

  character(len=*), parameter, private :: ModuleName='VENTFAC'

  public ventilation_3M, ventilation_1M_2M
contains

  subroutine ventilation_3M(ixy_inner, k, V, n0, lam, mu, params)

USE vernier_mod, ONLY: vernier_start, vernier_stop, vik
USE vernier_timing_mod, ONLY: LPROF_CASIM
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='VENTILATION_3M'

    integer, intent(in) :: ixy_inner
    integer, intent(in) :: k
    real(wp), intent(in) :: n0, lam, mu
    real(wp), intent(out) :: V  ! bulk ventilation factor
    type(hydro_params), intent(in) :: params

    real(wp) :: T1, T2
    real(wp) :: Sc ! Schmidt number
    real(wp) :: a_x, b_x, f_x

INTEGER(KIND=vik) :: vnr_handle_1

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (LPROF_CASIM) CALL vernier_start(vnr_handle_1, ModuleName//':'//RoutineName)

    a_x=params%a_x
    b_x=params%b_x
    f_x=params%f_x

    Sc=visair/Dv

    T1=vent_1*(mu+1.0)/lam
    T2=vent_2*Sc**(1.0/3.0)*(a_x*(rho0/rho(k,ixy_inner))**(0.5)*rho(k,ixy_inner)/visair)**(0.5)

!changing to rho on top give V units of 1/m2 which is correct to combine with melting, wet growth that have units (m2/s *V) to give 1/s for mass change rate
    V=2.0*pi*n0*rho(k,ixy_inner)*(T1+T2*GammaFunc(0.5*b_x+mu+2.5)/GammaFunc(1.0+mu) &
         *(1.0 + 0.5*f_x/lam)**(-(0.5*b_x + mu + 2.5))*lam**(-0.5*b_x-1.5))
!prf for single moment use capacitance of 0.5*sphere -i.e. ~plate or disk or aggregate
    if (l_kfsm) V=0.5*V
!prf

    ! for computational efficiency/accuracy changed from...
    ![' ']    !*(lam + 0.5*f_x)**(-(0.5*b_x + mu + 2.5))*lam**(1.+mu)) &

    IF (LPROF_CASIM) CALL vernier_stop(vnr_handle_1)

  end subroutine ventilation_3M

  subroutine ventilation_1M_2M(ixy_inner, k, V, n0, lam, mu, params)

USE vernier_mod, ONLY: vernier_start, vernier_stop, vik
USE vernier_timing_mod, ONLY: LPROF_CASIM
    USE parkind1, ONLY: jprb, jpim

    implicit none

    character(len=*), parameter :: RoutineName='VENTILATION_1M_2M'

    integer, intent(in) :: ixy_inner
    integer, intent(in) :: k
    real(wp), intent(in) :: n0, lam, mu
    real(wp), intent(out) :: V  ! bulk ventilation factor
    type(hydro_params), intent(in) :: params

    real(wp) :: T1, T2
    real(wp) :: Sc ! Schmidt number
    real(wp) :: a_x, b_x, f_x

INTEGER(KIND=vik) :: vnr_handle_2

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (LPROF_CASIM) CALL vernier_start(vnr_handle_2, ModuleName//':'//RoutineName)

    a_x=params%a_x
    b_x=params%b_x
    f_x=params%f_x

    Sc=visair/Dv

    T1=vent_1*(mu+1.0)/lam
    T2=vent_2*Sc**(1.0/3.0)*SQRT(a_x*SQRT(rho0*rho(k,ixy_inner))/visair)

    V=2.0*pi*n0*rho(k,ixy_inner)*(T1+T2*params%gam_0p5bx_mu_2p5/params%gam_1_mu &
         *(1.0 + 0.5*f_x/lam)**(-(0.5*b_x + mu + 2.5))*lam**(-0.5*b_x-1.5))
!prf for single moment use capacitance of 0.5*sphere -i.e. ~plate or disk or aggregate
    if (l_kfsm) V=0.5*V
!prf

    ! for computational efficiency/accuracy changed from...
    ![' ']    !*(lam + 0.5*f_x)**(-(0.5*b_x + mu + 2.5))*lam**(1.+mu)) &

    IF (LPROF_CASIM) CALL vernier_stop(vnr_handle_2)

  end subroutine ventilation_1M_2M
end module ventfac
