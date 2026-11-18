! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE init_ancillaries_coupling_mod

USE logging_mod, ONLY: log_fatal

IMPLICIT NONE

CONTAINS

#if defined(RIVER_CPL)
SUBROUTINE init_ancillaries_coupling(rivers, rivers_data)

USE jules_rivers_mod,         ONLY: nx_rivers, ny_rivers, l_riv_number,        &
                                    rivers_type, rivers_data_type
USE grid_utils_mod,           ONLY: grid_create, grid_info
USE oasis_rivers_mod,         ONLY: riv_number_file
USE oasis_rivers_control_mod, ONLY: n_rivers
USE missing_data_mod,         ONLY: rmdi
USE input_mod,                ONLY: dummy_grid => grid
USE fill_variables_from_file_mod, ONLY: fill_variables_from_file
!
!-------------------------------------------------------------------------------
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Fills river number ancillary if necessary
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! types
TYPE(rivers_type), INTENT(IN OUT) :: rivers
TYPE(rivers_data_type), INTENT(IN OUT), TARGET :: rivers_data

! Local variables
TYPE(grid_info) :: local_grid

INTEGER :: ERROR

CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'INIT_ANCILLARIES_COUPLING'

!------------------------------------------------------------------------------
! Read river number ancillary if necessary
IF (l_riv_number) THEN
  ! Temporarily copy saved grid (set for full model grid) to a local variable
  ! before overwriting to define the river number grid
  local_grid = grid_create(dummy_grid%is_1d,dummy_grid%dim_name,               &
                    dummy_grid%nx,dummy_grid%x_name,dummy_grid%nx,             &
                    dummy_grid%y_name,dummy_grid%ny)
  dummy_grid = grid_create( .FALSE., "", 0, 'longitude', nx_rivers,            &
                                             'latitude', ny_rivers)

  CALL fill_variables_from_file(TRIM(riv_number_file),                         &
                              [ 'rivers_outflow_number' ],[ 'river_number' ])

  n_rivers = MAXVAL(rivers%rivers_outflow_number(:,:))
  IF (n_rivers == rmdi) THEN
    CALL log_fatal("set_oasis_fields",                                         &
                   "river number ancillary contains missing data.")
  END IF

  ! Reset saved JULES grid from river to full land model
  dummy_grid = grid_create(local_grid%is_1d,local_grid%dim_name,               &
                           local_grid%nx,local_grid%x_name,local_grid%nx,      &
                           local_grid%y_name,local_grid%ny)
END IF

RETURN

END SUBROUTINE init_ancillaries_coupling
#else
SUBROUTINE init_ancillaries_coupling(rivers, rivers_data)
!
USE jules_rivers_mod, ONLY: rivers_type,rivers_data_type
!
!-------------------------------------------------------------------------------
IMPLICIT NONE

! types
TYPE(rivers_type), INTENT(IN OUT) :: rivers
TYPE(rivers_data_type), INTENT(IN OUT), TARGET :: rivers_data

! Local variables
CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'INIT_ANCILLARIES_COUPLING'

!-----------------------------------------------------------------------------
! Description: Stub routine for init_ancillaries_coupling
!              Run time controls should mean this routine is never
!              actually called.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

CALL log_fatal(RoutineName,                                                    &
                    "**ERROR**: init_ancillaries_coupling unavailable. " //    &
                    "Check RIVER_CPL fpp_defs key is set")

END SUBROUTINE init_ancillaries_coupling
#endif
END MODULE init_ancillaries_coupling_mod
