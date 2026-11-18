MODULE oasis_rivers_control_mod
! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! Please refer to the OASIS-MCT version 4.0 user guide for details on the implementation:
!  https://oasis.cerfacs.fr/wp-content/uploads/sites/114/2021/12/GLOBC_TR_oasis3mct_UserGuide_4.0_final_122021.pdf
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************
USE string_utils_mod, ONLY: to_string
USE logging_mod, ONLY: log_info, log_warn, log_fatal

IMPLICIT NONE

PRIVATE

! OASIS variables
INTEGER                   :: compid, il_part_id, il_part_id_1d
!
! Number of rivers in the river number ancillary file
INTEGER, PUBLIC           :: n_rivers
!
CHARACTER(LEN=5), PARAMETER :: comp_name='river'

! Public subroutines
PUBLIC :: oasis_init, oasis_partition, oasis_grid, oasis_variables,            &
          oasis_send, oasis_receive, oasis_finalise

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='OASIS_RIVERS_CONTROL_MOD'

CONTAINS

#if defined(RIVER_CPL)
SUBROUTINE oasis_init()

USE mod_oasis,                   ONLY: oasis_init_comp, oasis_get_localcomm
USE jules_vars_mod,              ONLY: mpi_local_comm
USE jules_model_environment_mod, ONLY: l_oasis_rivers

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises the oasis coupler and creates local MPI communicator
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Local variables
INTEGER :: ierror, ierror2

CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_INIT'

!------------------------------------------------------------------------------
! Coupling initialisation
CALL oasis_init_comp(compid, comp_name, ierror)

IF (ierror == 0) THEN
  ! Get the communicator allocated to this model component.
  CALL oasis_get_localcomm(mpi_local_comm, ierror2)
END IF

! Log messages need to be called after mpi_local_comm is initialized
IF (ierror /= 0) THEN
  CALL log_fatal(RoutineName,                                                  &
               "oasis_init_comp Error " //                                     &
               "(IOSTAT=" // TRIM(to_string(ierror)) // ")")
ELSE
  CALL log_info(RoutineName,"oasis_init_comp OK" )
END IF

IF (ierror2 /= 0) THEN
  CALL log_fatal(RoutineName,                                                  &
               "oasis_get_localcomm Error " //                                 &
               "(IOSTAT=" // TRIM(to_string(ierror2)) // ")")
ELSE
  CALL log_info(RoutineName,"oasis_get_localcomm OK" )
END IF

! Indicate that river coupling is taking place
l_oasis_rivers = .TRUE.

RETURN

END SUBROUTINE oasis_init
!
SUBROUTINE oasis_partition()

USE mod_oasis,        ONLY: oasis_def_partition
USE oasis_rivers_mod, ONLY: l_runoff_1d
USE jules_rivers_mod, ONLY: np_rivers
!
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines the oasis grid partition
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Local variables
INTEGER, DIMENSION(3) :: ig_paral
INTEGER :: ierror

CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_PARTITION'

!------------------------------------------------------------------------------
! Use serial partition, as rivers runs using only one processor
ig_paral(1) = 0
ig_paral(2) = 0

IF (l_runoff_1d) THEN
  ! Define the river outflow partition
  ig_paral(3) = n_rivers

  ! Coupling initialisation
  CALL oasis_def_partition(il_part_id_1d, ig_paral, ierror)

  IF (ierror /= 0) THEN
    CALL log_fatal(RoutineName,                                                &
                "oasis_def_partition 1d Error " //                             &
                "(IOSTAT=" // TRIM(to_string(ierror)) // ")")
  ELSE
    CALL log_info(RoutineName,"oasis_def_partition 1d OK" )
  END IF
END IF

ig_paral(3) = np_rivers

! Coupling initialisation
CALL oasis_def_partition(il_part_id, ig_paral, ierror)

IF (ierror /= 0) THEN
  CALL log_fatal(RoutineName,                                                  &
               "oasis_def_partition Error " //                                 &
               "(IOSTAT=" // TRIM(to_string(ierror)) // ")")
ELSE
  CALL log_info(RoutineName,"oasis_def_partition OK" )
END IF

RETURN

END SUBROUTINE oasis_partition
!
SUBROUTINE oasis_grid()
!
USE mod_oasis,        ONLY: oasis_start_grids_writing, oasis_write_grid,       &
                            oasis_write_corner, oasis_write_area,              &
                            oasis_write_mask, oasis_terminate_grids_writing
USE jules_rivers_mod, ONLY: np_rivers, rivers_dlat=>rivers_dy,                 &
                            rivers_dlon=>rivers_dx
USE jules_fields_mod, ONLY: rivers
!
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Writes to file the oasis grid
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Local variables
REAL    :: area(np_rivers,1),CornerLon(np_rivers,1,4),CornerLat(np_rivers,1,4)
INTEGER :: mask(np_rivers,1)
INTEGER :: i, ierror

CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_GRID'

!------------------------------------------------------------------------------
! Create grid file
CALL oasis_start_grids_writing(ierror)
!
! Calculate area, corners, mask: the river grid is regular
DO i = 1, np_rivers
  CornerLon(i,1,1) = rivers%rivers_lon_rp(i) - 0.5 * rivers_dlon
  CornerLon(i,1,2) = rivers%rivers_lon_rp(i) + 0.5 * rivers_dlon
  CornerLon(i,1,3) = CornerLon(i,1,2)
  CornerLon(i,1,4) = CornerLon(i,1,1)

  CornerLat(i,1,1) = rivers%rivers_lat_rp(i) - 0.5 * rivers_dlat
  CornerLat(i,1,2) = CornerLat(i,1,1)
  CornerLat(i,1,3) = rivers%rivers_lat_rp(i) + 0.5 * rivers_dlat
  CornerLat(i,1,4) = CornerLat(i,1,3)

  area(i,1) = rivers_dlat * rivers_dlon
  mask(i,1) = 1
END DO
!
CALL oasis_write_grid('rivt',np_rivers,1,                                      &
                      RESHAPE(rivers%rivers_lon_rp, [np_rivers,1 ]),           &
                      RESHAPE(rivers%rivers_lat_rp, [np_rivers,1 ]),il_part_id)
CALL oasis_write_corner('rivt',np_rivers,1,4,CornerLon,CornerLat,il_part_id)
CALL oasis_write_area('rivt',np_rivers,1,area,il_part_id)
CALL oasis_write_mask('rivt',np_rivers,1,mask,il_part_id)
!
! Close grid file
CALL oasis_terminate_grids_writing()

RETURN

END SUBROUTINE oasis_grid
!
SUBROUTINE oasis_variables()
!
USE mod_oasis,        ONLY: oasis_def_var, oasis_enddef, oasis_in, oasis_out,  &
                            oasis_real
USE oasis_rivers_mod, ONLY: np_receive, np_send, cpl_receive, cpl_send
!
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Declares all coupling variables
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Local variables
INTEGER, PARAMETER :: var_nodims(2) = [1,1 ]
INTEGER, PARAMETER :: var_actual_shape(2) = [1,1 ]
INTEGER :: i, ierror

CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_VARIABLES'

!------------------------------------------------------------------------------
! Declare oasis send fields
DO i = 1, np_send
  IF (TRIM(cpl_send(i)%field_name) == "rflow_outflow") THEN
    CALL oasis_def_var(cpl_send(i)%field_id, TRIM(cpl_send(i)%field_name),     &
                       il_part_id_1d, var_nodims, oasis_out, var_actual_shape, &
                       oasis_real, ierror)
  ELSE
    CALL oasis_def_var(cpl_send(i)%field_id, TRIM(cpl_send(i)%field_name),     &
                       il_part_id, var_nodims, oasis_out, var_actual_shape,    &
                       oasis_real, ierror)
  END IF
  IF (ierror /= 0) THEN
    CALL log_fatal(RoutineName,                                                &
                 "oasis_def_var Error " //                                     &
                 "(IOSTAT=" // TRIM(to_string(ierror)) // "): " //             &
                 "when defining variable " // TRIM(cpl_send(i)%field_name))
  END IF
  IF (cpl_send(i)%field_id == -1) THEN
    CALL log_fatal(RoutineName,                                                &
                 "oasis_def_var Error " //                                     &
                 "(IOSTAT=" // TRIM(to_string(ierror)) // "): " //             &
                 "variable " // TRIM(cpl_send(i)%field_name) //                &
                 " not found in the namcouple file")
  END IF
END DO

! Declare oasis receive fields
DO i = 1, np_receive
  CALL oasis_def_var(cpl_receive(i)%field_id, TRIM(cpl_receive(i)%field_name), &
                     il_part_id, var_nodims, oasis_in, var_actual_shape,       &
                     oasis_real, ierror)
  IF (ierror /= 0) THEN
    CALL log_fatal(RoutineName,                                                &
                 "oasis_def_var Error " //                                     &
                 "(IOSTAT=" // TRIM(to_string(ierror)) // "): " //             &
                 "when defining variable " // TRIM(cpl_receive(i)%field_name))
  END IF
  IF (cpl_receive(i)%field_id == -1) THEN
    CALL log_fatal(RoutineName,                                                &
                 "oasis_def_var Error " //                                     &
                 "(IOSTAT=" // TRIM(to_string(ierror)) // "): " //             &
                 "variable " // TRIM(cpl_send(i)%field_name) //                &
                 " not found in the namcouple file")
  END IF
END DO

! Finish oasis definitions
CALL oasis_enddef(ierror)

IF (ierror /= 0) THEN
  CALL log_fatal(RoutineName,                                                  &
               "oasis_enddef Error " //                                        &
               "(IOSTAT=" // TRIM(to_string(ierror)) // ")")
ELSE
  CALL log_info(RoutineName,"oasis_variables OK" )
END IF

RETURN

END SUBROUTINE oasis_variables
!
SUBROUTINE oasis_send(rivers, TIME)

USE mod_oasis,        ONLY: oasis_put, oasis_sent, oasis_loctrans,             &
                            oasis_torest, oasis_output, oasis_sentout,         &
                            oasis_torestout, oasis_waitgroup, oasis_ok
USE oasis_rivers_mod, ONLY: np_send, cpl_send
USE jules_rivers_mod, ONLY: rivers_type, l_riv_number, np_rivers
USE missing_data_mod, ONLY: imdi
!
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Sends fields to other models via oasis
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Arguments with intent(in)
TYPE(rivers_type), INTENT(IN) :: rivers
INTEGER, INTENT(IN) :: TIME

! Local variables
REAL, ALLOCATABLE :: rflow_outflow(:)
INTEGER :: i, ip, ierror, ERROR

CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_SEND'

!------------------------------------------------------------------------------
! Calculate the total runoff if necessary. This is done here so that other send
! variables do not have to wait for this to be calculated in the following loop
IF (l_riv_number) THEN
  ! Calculate river outflow for each river. This can be done like this here
  ! because the river executable runs in only one processor
  ALLOCATE(rflow_outflow(n_rivers), STAT = ERROR)
  rflow_outflow = 0.0
  DO ip = 1, np_rivers
    IF (rivers%rivers_outflow_number_rp(ip) > 0) THEN
      rflow_outflow(rivers%rivers_outflow_number_rp(ip)) =                     &
                        rflow_outflow(rivers%rivers_outflow_number_rp(ip)) +   &
                        rivers%riverout_rp(ip)
    END IF
  END DO
END IF

! Loop in all variables that are sent via coupling
DO i = 1, np_send
  SELECT CASE(TRIM(cpl_send(i)%field_name))
    ! send the right field
  CASE ('rflow_outflow')
    CALL oasis_put(cpl_send(i)%field_id,TIME,rflow_outflow,ierror)
  END SELECT

  IF (ierror /= oasis_sent .AND. ierror /= oasis_loctrans .AND.                &
      ierror /= oasis_torest .AND. ierror /= oasis_output  .AND.               &
      ierror /= oasis_sentout .AND. ierror /= oasis_torestout .AND.            &
      ierror /= oasis_waitgroup .AND. ierror /= oasis_ok) THEN
    CALL log_fatal(RoutineName,                                                &
                 "oasis_put Error " //                                         &
    "(IOSTAT=" // TRIM(to_string(ierror)) // "): " //                          &
    "when sending variable " // TRIM(cpl_send(i)%field_name))
  END IF
END DO

IF (l_riv_number) THEN
  DEALLOCATE(rflow_outflow)
END IF

RETURN

END SUBROUTINE oasis_send
!
SUBROUTINE oasis_receive(rivers, TIME)

USE mod_oasis,        ONLY: oasis_get, oasis_recvd, oasis_fromrest,            &
                            oasis_input, oasis_recvout, oasis_fromrestout,     &
                            oasis_ok
USE oasis_rivers_mod, ONLY: np_receive, cpl_receive
USE jules_rivers_mod, ONLY: rivers_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Receives fields from other models via oasis
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Arguments with intent(in)
TYPE(rivers_type), INTENT(IN) :: rivers
INTEGER, INTENT(IN) :: TIME

! Local variables
INTEGER :: i, ierror

CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_RECEIVE'

!------------------------------------------------------------------------------
! Loop in all variables that are received via coupling
DO i = 1, np_receive
  SELECT CASE(TRIM(cpl_receive(i)%field_name))
    ! receive the right field
  CASE ('sub_surf_roff_rp')
    CALL oasis_get(cpl_receive(i)%field_id,TIME,rivers%sub_surf_roff_rp,ierror)
  CASE ('surf_roff_rp')
    CALL oasis_get(cpl_receive(i)%field_id,TIME,rivers%surf_roff_rp,ierror)
  END SELECT

  IF (ierror /= oasis_recvd .AND. ierror /= oasis_fromrest .AND.               &
      ierror /= oasis_input .AND. ierror /= oasis_recvout  .AND.               &
      ierror /= oasis_fromrestout .AND. ierror /= oasis_ok ) THEN
    CALL log_fatal(RoutineName,                                                &
                 "oasis_get Error " //                                         &
    "(IOSTAT=" // TRIM(to_string(ierror)) // "): " //                          &
    "when receiving variable " // TRIM(cpl_receive(i)%field_name))
  END IF
END DO

RETURN

END SUBROUTINE oasis_receive
!
SUBROUTINE oasis_finalise()
!
USE mod_oasis,        ONLY: oasis_terminate
!
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Finalises the oasis coupler
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Local variables
INTEGER :: ierror

CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_FINALISE'

!------------------------------------------------------------------------------
! Coupling finalisation
CALL oasis_terminate(ierror)

IF (ierror /= 0) THEN
  CALL log_fatal(RoutineName,                                                  &
               "oasis_terminate Error " //                                     &
               "(IOSTAT=" // TRIM(to_string(ierror)) // ")")
ELSE
  CALL log_info(RoutineName,"oasis_terminate OK" )
END IF

RETURN

END SUBROUTINE oasis_finalise
#else
SUBROUTINE oasis_init()

USE jules_model_environment_mod, ONLY: l_oasis_rivers

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description: Stub routine for oasis_init.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Indicate that river coupling is not taking place
l_oasis_rivers = .FALSE.

END SUBROUTINE oasis_init

SUBROUTINE oasis_partition()
!
IMPLICIT NONE

! Local variables
CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_PARTITION_STUB'

!-----------------------------------------------------------------------------
! Description: Stub routine for oasis_partition.
!              Run time controls should mean this routine is never
!              actually called.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

CALL log_fatal(RoutineName, "**ERROR**: oasis_partition unavailable. " //      &
                         "Check RIVER_CPL fpp_defs key is set")

END SUBROUTINE oasis_partition

SUBROUTINE oasis_grid()
!
IMPLICIT NONE

! Local variables
CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_GRID_STUB'

!-----------------------------------------------------------------------------
! Description: Stub routine for oasis_grid.
!              Run time controls should mean this routine is never
!              actually called.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

CALL log_fatal(RoutineName, "**ERROR**: oasis_grid unavailable. " //           &
                         "Check RIVER_CPL fpp_defs key is set")

END SUBROUTINE oasis_grid

SUBROUTINE set_oasis_fields()
!
IMPLICIT NONE

! Local variables
CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'SET_OASIS_FIELDS'

!-----------------------------------------------------------------------------
! Description: Stub routine for set_oasis_fields.
!              Run time controls should mean this routine is never
!              actually called.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

CALL log_fatal(RoutineName, "**ERROR**: set_oasis_fields unavailable. " //     &
                         "Check RIVER_CPL fpp_defs key is set")

END SUBROUTINE set_oasis_fields

SUBROUTINE oasis_variables()
!
IMPLICIT NONE

! Local variables
CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_VARIABLES'

!-----------------------------------------------------------------------------
! Description: Stub routine for oasis_variables.
!              Run time controls should mean this routine is never
!              actually called.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

CALL log_fatal(RoutineName, "**ERROR**: oasis_variables unavailable. " //      &
                         "Check RIVER_CPL fpp_defs key is set")


END SUBROUTINE oasis_variables

SUBROUTINE oasis_send(rivers,TIME)
!
USE jules_rivers_mod, ONLY: rivers_type,rivers_data_type
!
IMPLICIT NONE

! Arguments with intent(in)
TYPE(rivers_type), INTENT(IN OUT) :: rivers
INTEGER, INTENT(IN) :: TIME

! Local variables
CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_SEND'

!-----------------------------------------------------------------------------
! Description: Stub routine for oasis_send.
!              Run time controls should mean this routine is never
!              actually called.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

CALL log_fatal(RoutineName, "**ERROR**: oasis_send unavailable. " //           &
                         "Check RIVER_CPL fpp_defs key is set")

END SUBROUTINE oasis_send

SUBROUTINE oasis_receive(rivers,TIME)
!
USE jules_rivers_mod, ONLY: rivers_type,rivers_data_type
!
IMPLICIT NONE

! Arguments with intent(in)
TYPE(rivers_type), INTENT(IN OUT) :: rivers
INTEGER, INTENT(IN) :: TIME

! Local variables
CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_RECEIVE'

!-----------------------------------------------------------------------------
! Description: Stub routine for oasis_receive.
!              Run time controls should mean this routine is never
!              actually called.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

CALL log_fatal(RoutineName, "**ERROR**: oasis_receive unavailable. " //        &
                         "Check RIVER_CPL fpp_defs key is set")

END SUBROUTINE oasis_receive

SUBROUTINE oasis_finalise()
!
IMPLICIT NONE

! Local variables
CHARACTER(LEN=*), PARAMETER  :: RoutineName = 'OASIS_FINALISE'

!-----------------------------------------------------------------------------
! Description: Stub routine for oasis_finalise.
!              Run time controls should mean this routine is never
!              actually called.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

CALL log_fatal(RoutineName, "**ERROR**: oasis_finalise unavailable. " //       &
                         "Check RIVER_CPL fpp_defs key is set")

END SUBROUTINE oasis_finalise
#endif

END MODULE oasis_rivers_control_mod
