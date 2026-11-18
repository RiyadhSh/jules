!******************************COPYRIGHT****************************************
! (c) British Antarctic Survey and University of Bristol, 2023.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT****************************************

MODULE jules_water_tracers_mod

!-------------------------------------------------------------------------------
! Description:
!   This module contains water tracer settings and common functions
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in Surface
!
! Code Description:
!   Language: Fortran 90.
!-------------------------------------------------------------------------------

USE missing_data_mod, ONLY: imdi, rmdi
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! Switch to control water tracers and set the number of water tracers in JULES

! Note, water tracers are not available in standalone JULES.
! Also, the number of water tracers used in JULES can be a subset of the
! number used in UM.
! The first water tracer is always the normal water tracer.

#if defined(UM_JULES)
! **In the future, these will be passed in to JULES from the UM.
LOGICAL, PARAMETER :: l_wtrac_jls = .FALSE.
INTEGER, PARAMETER :: n_wtrac_jls = 1
#else
! Water tracers are just dummy fields in standalone JULES
LOGICAL, PARAMETER :: l_wtrac_jls = .FALSE.
INTEGER, PARAMETER :: n_wtrac_jls = 1
#endif

! Logical to control the water tracer calculations on the implicit scheme
! (Water tracers only need to be updated on the final dynamics loop in
!  UMJULES)
! **In the future, this logical will be passed into JULES from the UM in
!   the call to surf_couple_implicit
LOGICAL :: l_wtrac_imp_jls = .FALSE.

! Number of evaporative sources modelled:
!   snow, canopy intercepted water, soil evapotranspiration (bare soil
!   evap and transpiration from vegetative tiles),  and open water
!   (i.e. lakes).

INTEGER, PARAMETER :: n_evap_srce = 4

INTEGER, PARAMETER :: evap_srce_snow = 1
INTEGER, PARAMETER :: evap_srce_cano = 2
INTEGER, PARAMETER :: evap_srce_soil = 3
INTEGER, PARAMETER :: evap_srce_lake = 4

! Fields that will be set in a namelist in the future:

! Standard ratio for water tracers
REAL(KIND=real_jlslsm), PARAMETER :: standard_ratio_wtrac(n_wtrac_jls) = 1.0

! Water tracer fixed ratio for lakes
REAL(KIND=real_jlslsm), PARAMETER :: wtrac_lake_fixed_ratio = 1.0

! Minimum water value used to calculate water tracer ratios
REAL(KIND=real_jlslsm), PARAMETER :: min_q_ratio = 1.0e-18

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
                            ModuleName='JULES_WATER_TRACERS_MOD'

CONTAINS

!----------------------------------------------------------------------------

REAL FUNCTION wtrac_calc_ratio_fn_jules(i_wt, q_wtrac, q)

!-----------------------------------------------------------------------------
! Description:
!
!  Function to calculate the ratio of water tracer to normal water.
!  If there is no water, set ratio to Standard value.
!  This function is for 'real_jlslsm' precision inputs/outputs (i.e. used
!  in JULES).
!
!-----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER, INTENT(IN) :: i_wt                     ! Water tracer number

REAL(KIND=real_jlslsm), INTENT(IN) :: q_wtrac   ! Water tracer
REAL(KIND=real_jlslsm), INTENT(IN) :: q         ! Normal water

REAL(KIND=real_jlslsm) :: q_ratio               ! Output ratio
REAL(KIND=real_jlslsm) :: q_ratio_std           ! Standard ratio

q_ratio_std = standard_ratio_wtrac(i_wt)

IF (ABS(q) > min_q_ratio) THEN
  q_ratio = q_wtrac / q
ELSE
  q_ratio = q_ratio_std
END IF

! Yet to be decided if it's best to block negative ratios of not.
!IF (q_ratio < 0.0) q_ratio = q_ratio_std

wtrac_calc_ratio_fn_jules = q_ratio

END FUNCTION wtrac_calc_ratio_fn_jules

END MODULE jules_water_tracers_mod
