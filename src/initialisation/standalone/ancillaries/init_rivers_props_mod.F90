#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE init_rivers_props_mod

USE jules_rivers_mod, ONLY: inland_drainage, river_mouth

USE logging_mod, ONLY: log_fatal, log_info, log_warn

USE missing_data_mod, ONLY: imdi, rmdi

USE string_utils_mod, ONLY: to_string

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!------------------------------------------------------------------------------
! Scalar parameters.
!------------------------------------------------------------------------------
INTEGER, PARAMETER :: flow_dir_sea = 0
    ! The value of flow direction that is used to indicate
    ! a sea point. This value is used to set the flow
    ! direction variable in this code, but it need not be
    ! the value used to indicate sea in the input ancillary
    ! file. This value must not be any of 1:8, river_mouth,
    ! inland_drainage, or land_off_domain.

INTEGER, PARAMETER :: land_off_domain = -1
    ! The value of flow direction that is used to indicate
    ! a land point (according to the river ancillary) that is
    ! outside of the valid river domain, i.e. it is on the
    ! river grid but not in the area covered by land, hence
    ! it cannot become a river point.
    ! This value must not be any of 1:8, river_mouth,
    ! inland_drainage, or flow_dir_sea.
    ! Note that although this value is used to set the flow
    ! direction for points n the edge of the domain, those
    ! points are not currently accessed again - but in
    ! future it might be useful to know whether they flow to
    ! sea or another river point.

REAL(KIND=real_jlslsm), PARAMETER :: frac_toler = 1.0e-3
  ! Tolerance (fraction of a gridbox) used when testing whether model points
  ! lie on a regular grid, if grids have the same resolution, and when setting
  ! some bounding values.

!------------------------------------------------------------------------------
! Array parameters.
!------------------------------------------------------------------------------
INTEGER, PARAMETER :: flow_dir_river(10) =                                     &
  [ 1, 2, 3, 4, 5, 6, 7, 8, river_mouth, inland_drainage ]
    ! The values of flow direction that are acceptable for a river point.
    ! Effectively these also the only values that indicate land points in the
    ! flow direction field. The first 8 places represent "compass directions"
    ! for the 8 points immediately neighbouring a given point.

INTEGER, ALLOCATABLE :: grid_riv_pt_number(:,:)
  ! Map full river grid to river points. For each location on the grid, this
  ! is the river point number.

INTEGER, ALLOCATABLE :: rivers_dir_int(:,:)
  ! Integer version of rivers_dir (the river routing direction index).
  ! We use the REAL version for i.o (JULES cannot read an integer ancillary
  ! field) but thereafter the integer version is easier to use.

PRIVATE  !  private scope by default
PUBLIC init_rivers_props

CONTAINS

!##############################################################################

SUBROUTINE init_rivers_props(rivers, rivers_data)

USE ancil_namelist_mod, ONLY: check_namelist_values, initialise_namelist_vars, &
                              set_constant_vars, set_file_vars

USE fill_variables_from_file_mod, ONLY: fill_variables_from_file

USE grid_utils_mod, ONLY: grid_create, grid_info

USE init_ancillaries_coupling_mod, ONLY: init_ancillaries_coupling

USE input_mod, ONLY:                                                           &
    od_grid => grid,                                                           &
    dummy_grid => grid,                                                        &
    use_subgrid

USE init_overbank_props_mod, ONLY: init_overbank_props

USE io_constants, ONLY: max_sdf_name_len, max_file_name_len, namelist_unit

USE jules_model_environment_mod, ONLY: l_oasis_rivers

USE jules_rivers_mod, ONLY:                                                    &
  i_river_vn, land_dx, land_dy, l_rivers, nx_land_grid,                        &
  nx_rivers_in=>nx_rivers, ny_land_grid, ny_rivers_in=>ny_rivers,              &
  rivers_length, rivers_regrid, rivers_rfm, rivers_trip, x1_land_grid,         &
  y1_land_grid, l_riv_overbank,                                                &
  ! types
  rivers_type, rivers_data_type

USE model_grid_mod, ONLY: global_land_pts, l_coord_latlon

USE model_interface_mod, ONLY: identifier_len

USE overbank_inundation_mod, ONLY:                                             &
  init_rosgen_vars, overbank_model, overbank_simple_rosgen

USE parallel_mod, ONLY: is_master_task

USE templating_mod, ONLY: tpl_has_var_name

IMPLICIT NONE

!------------------------------------------------------------------------------
! Description:
!    Top-level routine that controls reading of river routing grid and river
!    and properties, and allocation and setting of river point values.
!
!    When running the river executable in coupled mode, no regridding between
!    the land and river grids is necessary, so we can skip the calculation of
!    some variables that will not be used.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Subroutine arguments.
!------------------------------------------------------------------------------
! Types
TYPE(rivers_type), INTENT(IN OUT) :: rivers
TYPE(rivers_data_type), INTENT(IN OUT), TARGET :: rivers_data

!------------------------------------------------------------------------------
! Local parameters.
!------------------------------------------------------------------------------
INTEGER, PARAMETER :: nvars_max = 6
  ! The maximum possible number of routing variables.

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'INIT_RIVERS_PROPS'

!------------------------------------------------------------------------------
! Local variables appearing in the namelist.
!------------------------------------------------------------------------------
! Scalars in namelist.
INTEGER :: nvars
    ! The number of variables in this section

INTEGER :: nx_rivers, ny_rivers
    ! The sizes of the x and y dimensions of the river routing input grid.

LOGICAL ::                                                                     &
  l_use_area = .FALSE.
    ! Switch (used with RFM) to use a drainage area ancillary field to identify
    ! river points.
    ! T means use ancillary data for area.
    ! F means no ancillary data provided. Area will be set above threshold so
    ! that all points will be identified as rivers.

CHARACTER(LEN=max_file_name_len) :: FILE
    ! The name of the file (or variable name template) to
    ! use for variables that need to be filled from file.
    ! Co-ordinate variables are also read from this file,
    ! assuming file does not include variable-name
    ! templating (i.e. assuming all ancillary variables will
    ! come from the same file).

CHARACTER(LEN=max_file_name_len) :: coordinate_file
    ! The name of the file given for reading coordinate
    ! values from. This is only used if the ancillary file
    ! name includes templating (i.e. ancillary variables will
    ! come from multiple files), in which case this variable
    ! is used to provide clarity as to where the coordinates
    ! are read from.

CHARACTER(LEN=max_sdf_name_len) :: x_dim_name, y_dim_name
    ! The names of the x and y dimensions.

! Arrays in namelist.
REAL(KIND=real_jlslsm) :: const_val(nvars_max)
    ! The constant value to use for each variable if
    ! use_file = F for that variable

LOGICAL :: use_file(nvars_max)
    !   T - the variable uses the file
    !   F - the variable is set using a constant value
LOGICAL ::                                                                     &
  use_var(nvars_max)
    ! Flag indicating if a variable is required.

CHARACTER(LEN=max_sdf_name_len) :: tpl_name(nvars_max)
    ! The name to substitute in a template for each variable
CHARACTER(LEN=identifier_len) :: var(nvars_max)
    ! The variable identifiers of the variables
CHARACTER(LEN=max_sdf_name_len) :: var_name(nvars_max)
    ! The name of each variable in the file

!------------------------------------------------------------------------------
! Definition of the namelist. This contains local variables and some from
! jules_rivers_mod.
!------------------------------------------------------------------------------
NAMELIST  /jules_rivers_props/                                                 &
   x_dim_name, y_dim_name, nx_rivers, ny_rivers, FILE, coordinate_file,        &
   nvars, var, use_file, var_name, tpl_name, const_val,                        &
   land_dx, land_dy, l_use_area, nx_land_grid, ny_land_grid, rivers_length,    &
   rivers_regrid, x1_land_grid, y1_land_grid

!------------------------------------------------------------------------------
! Local scalar variables (not in namelist).
!------------------------------------------------------------------------------
INTEGER :: ERROR, error_sum    ! Error indicators

INTEGER :: nvars_const         ! The number of variables to be set using
                               ! constant values.
INTEGER :: nvars_file          ! The number of variables that will be set
                               ! from the given file (template?)
INTEGER :: nvars_required      ! The number of routing variables that are
                               ! required in this configuration

LOGICAL ::                                                                     &
  l_calc_rivers_length,                                                        &
    ! Flag indicating if rivers_length will be calculated.
  l_have_template,                                                             &
    ! Flag indicating if variable-name templating is being used.
  use_sub_local
    ! A local copy of use_subgrid.

!------------------------------------------------------------------------------
! Local array variables (not in namelist).
!------------------------------------------------------------------------------
CHARACTER(LEN=max_file_name_len) :: file_name(nvars_max)
    ! The name of the file to use for each variable.
CHARACTER(LEN=max_file_name_len) :: file_name_coords
    ! The name of the file to be used to read coordinate values from.
CHARACTER(LEN=identifier_len) :: required_vars(nvars_max)
    ! The variable identifiers of the required variables

!------------------------------------------------------------------------------
! Types
!------------------------------------------------------------------------------
TYPE(grid_info) :: local_grid

!end of header
!------------------------------------------------------------------------------

! If rivers are not on, there is nothing to do.
IF ( .NOT. l_rivers ) RETURN

! If this run does not have any land points we cannot do river routing.
IF ( global_land_pts == 0 ) THEN
  CALL log_fatal( RoutineName,                                                 &
                 "River routing not appropriate for runs with no land "     // &
                 "points." )
END IF

!------------------------------------------------------------------------------
! Most tasks need only be done on the master task, which is where the river
! routing code will run. If the JULES_RIVERS_PROPS namelist were to contain
! values that are needed by others tasks, all tasks would have to read it.
!------------------------------------------------------------------------------
IF ( is_master_task() ) THEN

  !----------------------------------------------------------------------------
  ! Get ready for reading the river routing properties namelist.
  !----------------------------------------------------------------------------
  ! Initialisation of some common variables related to namelists.
  CALL initialise_namelist_vars( nvars, nvars_const, nvars_file,               &
                                 nvars_required, l_have_template,              &
                                 FILE, const_val, use_file, use_var,           &
                                 file_name, tpl_name, var, var_name )

  ! Initialisation of further variables.
  coordinate_file = '' ! Empty file name.
  l_calc_rivers_length = .FALSE.  !  value will not be calculated

  !----------------------------------------------------------------------------
  ! Read river routing properties namelist
  !----------------------------------------------------------------------------
  CALL log_info( RoutineName, "Reading JULES_RIVERS_PROPS namelist..." )

  READ(namelist_unit, NML = jules_rivers_props, IOSTAT = ERROR)
  IF ( ERROR /= 0 ) THEN
    CALL log_fatal( RoutineName,                                               &
                    "Error reading namelist JULES_RIVERS_PROPS "            // &
                    "(IOSTAT=" // TRIM(to_string(ERROR)) // ")" )
  END IF

  !----------------------------------------------------------------------------
  ! Pass the values of nx_rivers and ny_rivers that were read from
  ! jules_rivers_props to the variables in jules_rivers_mod.
  ! Note that these are only needed by the master task.
  !----------------------------------------------------------------------------
  nx_rivers_in = nx_rivers
  ny_rivers_in = ny_rivers

  !----------------------------------------------------------------------------
  ! Check values from the namelist.
  !----------------------------------------------------------------------------
  CALL check_jules_rivers_props( nx_rivers, ny_rivers, coordinate_file, FILE,  &
                                 l_calc_rivers_length )

  !----------------------------------------------------------------------------
  ! Set up the list of required ancillary variables.
  !----------------------------------------------------------------------------
  ! Add variables that depend only on the chosen routing model.
  SELECT CASE ( i_river_vn )

  CASE ( rivers_trip )

    nvars_required = 2
    required_vars(1:nvars_required) = [ 'direction', 'sequence ' ]

  CASE ( rivers_rfm )

    nvars_required = 1
    required_vars(1:nvars_required) = [ 'direction' ]

  CASE DEFAULT
    CALL log_fatal( RoutineName,                                               &
                    "Do not recognise i_river_vn = '"                       // &
                    TRIM(to_string(i_river_vn)) // "'" )
  END SELECT

  ! Add required variables that depend on the coordinates used.
  ! If the river grid is not a lat-lon grid, we need lat and lon variables.
  ! These are needed for output files, and for calculating gridbox areas
  ! under some configurations.
  IF ( .NOT. l_coord_latlon ) THEN
    nvars_required = nvars_required + 2
    required_vars(nvars_required-1:nvars_required) =                           &
      [ 'latitude_2d ', 'longitude_2d' ]
  END IF

  ! Add variables that are required by some configurations.
  IF ( i_river_vn == rivers_rfm .AND. l_use_area ) THEN
    ! Add upstream drainage area.
    nvars_required = nvars_required + 1
    required_vars(nvars_required) = 'area'
  END IF

  !----------------------------------------------------------------------------
  ! Check that variables in the namelist have been provided and set up other
  ! variables.
  !----------------------------------------------------------------------------
  CALL check_namelist_values( nvars, nvars_max, nvars_required,                &
                              RoutineName, FILE, const_val,                    &
                              use_file, required_vars, tpl_name, var,          &
                              nvars_const, nvars_file, l_have_template,        &
                              use_var, file_name, var_name )

  !----------------------------------------------------------------------------
  ! Set up the regular routing input grid.
  !----------------------------------------------------------------------------
  use_sub_local = use_subgrid
  use_subgrid   = .FALSE.

  ! Temporarily copy saved grid (set for full model grid) to a local variable
  ! before overwriting to define the river routing grid
  local_grid = grid_create( dummy_grid%is_1d, dummy_grid%dim_name,             &
                            dummy_grid%nx, dummy_grid%x_name, dummy_grid%nx,   &
                            dummy_grid%y_name, dummy_grid%ny )

  ! Get the name of the file to be used to read coordinates. This depends on
  ! whether variable-name templating is used. We have already checked that a
  ! file name was provided.

  ! Establish if templating is specified.
  IF ( tpl_has_var_name(FILE) ) THEN
    ! We use a separate variable to get the file name.
    ! We have previously checked that a file name was provided.
    file_name_coords = coordinate_file
  ELSE
    ! Use the only ancillary file.
    file_name_coords = FILE
    ! If another file was indicated, clarify that it will not be used.
    IF ( LEN_TRIM(coordinate_file) > 0 ) THEN
      CALL log_info( RoutineName,                                              &
                     "No templating; coordinate_file will be ignored.")
    END IF
  END IF

  ! Read y dimension from file.
  od_grid = grid_create( grid_is_1d = .TRUE.,                                  &
                         dim_name = TRIM(y_dim_name), npoints = ny_rivers )
  ALLOCATE(rivers_data%rivers_ygrid(ny_rivers), STAT = ERROR)
  rivers%rivers_ygrid => rivers_data%rivers_ygrid
  CALL fill_variables_from_file(file_name_coords,                              &
                                [ 'rivers_ygrid' ],[ TRIM(y_dim_name) ])

  ! Read x dimension from file.
  od_grid = grid_create( grid_is_1d = .TRUE.,                                  &
                         dim_name = TRIM(x_dim_name), npoints = nx_rivers)
  ALLOCATE(rivers_data%rivers_xgrid(nx_rivers), STAT = ERROR)
  rivers%rivers_xgrid => rivers_data%rivers_xgrid
  CALL fill_variables_from_file(file_name_coords,                              &
                                [ 'rivers_xgrid' ],[ TRIM(x_dim_name) ])

  !----------------------------------------------------------------------------
  ! Build the routing grid object from the namelist values
  !----------------------------------------------------------------------------
  dummy_grid = grid_create( .FALSE., "", 0, x_dim_name, nx_rivers, y_dim_name, &
                            ny_rivers)

  !----------------------------------------------------------------------------
  ! Allocate gridded river ancillary variables.
  !----------------------------------------------------------------------------
  CALL allocate_river_vars_grid( global_land_pts, nx_rivers, ny_rivers, rivers,&
                                 rivers_data )

  !----------------------------------------------------------------------------
  ! Set variables that are using a constant value.
  !----------------------------------------------------------------------------
  CALL set_constant_vars( nvars, nvars_const, nvars_max, const_val,            &
                          use_file, use_var, var )

  !----------------------------------------------------------------------------
  ! Set variables that are to be read from a file.
  !----------------------------------------------------------------------------
  CALL set_file_vars( nvars, nvars_file, nvars_max, l_have_template,           &
                      RoutineName, use_file, use_var, file_name, var, var_name)

END IF  !  is_master_task

!-----------------------------------------------------------------------------
! Read overbank and rivers coupling properties namelists if required - so that
! any required ancillary fields are read from file before translating to the
! river routing points.
!-----------------------------------------------------------------------------
IF ( l_riv_overbank ) THEN
  CALL init_overbank_props()
END IF

IF ( l_oasis_rivers ) THEN
  CALL init_ancillaries_coupling( rivers, rivers_data )
END IF

!------------------------------------------------------------------------------
! Process the ancillary fields to identify river points and populate fields
! on river points.
!------------------------------------------------------------------------------
CALL process_rivers_data( l_calc_rivers_length, l_use_area,                    &
                          rivers, rivers_data )

IF ( is_master_task() ) THEN

  !----------------------------------------------------------------------------
  ! Calculate derived variables for inundation.
  !----------------------------------------------------------------------------
  IF ( l_riv_overbank .AND. overbank_model == overbank_simple_rosgen ) THEN
    CALL init_rosgen_vars( rivers )
  END IF

  !----------------------------------------------------------------------------
  ! Reset saved JULES grid from river to full land model
  !----------------------------------------------------------------------------
  dummy_grid = grid_create(local_grid%is_1d,local_grid%dim_name,               &
                           local_grid%nx,local_grid%x_name,local_grid%nx,      &
                           local_grid%y_name,local_grid%ny)
  use_subgrid = use_sub_local

END IF  !  is_master_task

RETURN
END SUBROUTINE init_rivers_props

!##############################################################################

SUBROUTINE check_jules_rivers_props( nx_rivers, ny_rivers, coordinate_file,    &
                                     FILE, l_calc_rivers_length )

USE jules_rivers_mod, ONLY: i_river_vn, land_dx, land_dy, l_riv_overbank,      &
                            nx_land_grid, ny_land_grid, x1_land_grid,          &
                            y1_land_grid, rivers_length,                       &
                            rivers_regrid, rivers_rfm

USE model_grid_mod, ONLY: l_coord_latlon

USE overbank_inundation_mod, ONLY: overbank_simple, overbank_simple_rosgen,    &
                                   overbank_model

USE templating_mod, ONLY: tpl_has_var_name

IMPLICIT NONE

!------------------------------------------------------------------------------
! Scalar arguments with INTENT(IN)
!------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  nx_rivers, ny_rivers
    ! The sizes of the x and y dimensions of the river routing input grid.

CHARACTER(LEN=*), INTENT(IN) ::                                                &
  coordinate_file,                                                             &
    ! The name of the file given for reading coordinate
    ! values from.
  FILE
    ! The name of the file to be used to read data from.

!------------------------------------------------------------------------------
! Scalar arguments with INTENT(IN OUT)
!------------------------------------------------------------------------------
LOGICAL, INTENT(IN OUT) ::                                                     &
  l_calc_rivers_length
    ! Flag indicating if rivers_length will be calculated.

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CHECK_JULES_RIVERS_PROPS'

!end of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Check that a file name was provided.
!------------------------------------------------------------------------------
IF ( LEN_TRIM(FILE) == 0 ) THEN
  CALL log_fatal( RoutineName, "No name given for file." )
END IF

! If templating is specified, check we have a coordinate_file.
IF ( tpl_has_var_name(FILE) .AND. LEN_TRIM(coordinate_file) == 0 ) THEN
  CALL log_fatal( RoutineName, "No name given for coordinate_file." )
END IF

!------------------------------------------------------------------------------
! Regridding is only supported for (regular) lat-lon grids.
! The regularity is checked later.
!------------------------------------------------------------------------------
IF ( rivers_regrid .AND. .NOT. l_coord_latlon ) THEN
  CALL log_fatal( RoutineName,                                                 &
                 "Regridding can only be used with lat-lon grids." )
END IF

!------------------------------------------------------------------------------
! Insist on a river grid that is at least 2 x 2. This is required for the
! calculations of rivers_dx and rivers_dy.
!------------------------------------------------------------------------------
IF ( nx_rivers < 2  .OR. ny_rivers < 2 ) THEN
  CALL log_fatal( RoutineName,                                                 &
                  "Routing grid must be at least 2x2. Dimensions provided: " //&
                  "nx_rivers = " // TRIM(to_string(nx_rivers))               //&
                  "ny_rivers = " // TRIM(to_string(ny_rivers)) )
END IF

!------------------------------------------------------------------------------
! land_dy and land_dx are required, and should be >0.
!------------------------------------------------------------------------------
IF ( ABS( land_dy - rmdi ) < EPSILON(1.0) .OR.                                 &
     ABS( land_dx - rmdi ) < EPSILON(1.0) ) THEN
  CALL log_fatal( RoutineName,                                                 &
                  "land_dy and land_dx must be provided." )
ELSE IF ( land_dx <= 0.0 .OR. land_dy <= 0.0 ) THEN
  CALL log_fatal( RoutineName,                                                 &
                  "land_dy and land_dx must be >0." )
END IF

!------------------------------------------------------------------------------
! Conditions for nx_land_grid, ny_land_grid, x1_land_grid and y1_land_grid.
! Note that we later test that nx_land_grid and ny_land_grid are large enough.
!------------------------------------------------------------------------------
IF ( nx_land_grid == imdi .OR. ny_land_grid == imdi .OR.                       &
     ABS( y1_land_grid - rmdi ) < EPSILON(1.0)  .OR.                           &
     ABS( x1_land_grid - rmdi ) < EPSILON(1.0) ) THEN
  CALL log_fatal( RoutineName,                                                 &
                  "nx_land_grid, ny_land_grid, x1_land_grid and "           // &
                  "y1_land_grid must be provided for this configuration." )
ELSE IF ( nx_land_grid <= 0 .OR. ny_land_grid <= 0 ) THEN
  CALL log_fatal( RoutineName,                                                 &
                  "nx_land_grid and ny_land_grid must be >0.")
END IF

!------------------------------------------------------------------------------
! Conditions for rivers_length.
!------------------------------------------------------------------------------
IF ( .NOT. l_coord_latlon ) THEN
  ! We don't have a latitude-longitude grid, and so must have rivers_length>0
  ! because that is used to calculate box areas.
  IF ( ABS( rivers_length - rmdi ) < EPSILON(1.0) .OR.                         &
       rivers_length <= 0.0 ) THEN
    CALL log_fatal( RoutineName,                                               &
                    "rivers_length > 0 must be provided for this "          // &
                    "configuration." )
  END IF
END IF

! Note that this check requires that we have read the jules_overbank namelist.
IF ( i_river_vn == rivers_rfm .OR.                                             &
     ( l_riv_overbank .AND. ( overbank_model == overbank_simple .OR.           &
                              overbank_model == overbank_simple_rosgen ) )     &
   ) THEN
  ! These configurations require rivers_length. A value <=0 will later trigger
  ! calculation of rivers_length.
  IF ( ABS( rivers_length - rmdi ) < EPSILON(1.0) ) THEN
    CALL log_fatal( RoutineName,                                               &
                    "rivers_length must be provided for this configuration." )
  ELSE IF ( rivers_length <= 0.0 ) THEN
    ! A non-missing value has been provided but it is <= 0, hence a new value
    ! will be calculated - this is so as to mimic the behaviour of earlier
    ! versions of the code. Note that we will only get here for grids that are
    ! using lat and lon coords, because non-latlon grids were tested earlier
    ! to ensure rivers_length > 0.
    ! Set a switch to trigger later calculation of rivers_length.
    l_calc_rivers_length = .TRUE.
    CALL log_info( RoutineName,                                                &
                   "rivers_length <= 0: value will be calculated from grid." )
  END IF
END IF

RETURN
END SUBROUTINE check_jules_rivers_props

!##############################################################################

SUBROUTINE allocate_river_vars_grid( global_land_pts, nx_rivers, ny_rivers,    &
                                     rivers, rivers_data )

!------------------------------------------------------------------------------
! Description:
!   Allocate gridded river ancillary variables, initialise, and associate
!   pointers.
!------------------------------------------------------------------------------

USE jules_rivers_mod, ONLY: i_river_vn, l_riv_number, rivers_data_type,        &
                            rivers_rfm, rivers_trip, rivers_type

IMPLICIT NONE

!------------------------------------------------------------------------------
! Scalar arguments with INTENT(IN)
!------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  global_land_pts,                                                             &
    ! Number of land points.
  nx_rivers, ny_rivers
    ! Grid sizes.

!------------------------------------------------------------------------------
! Arguments with INTENT(IN OUT)
!------------------------------------------------------------------------------
TYPE(rivers_type), INTENT(IN OUT) :: rivers
TYPE(rivers_data_type), INTENT(IN OUT), TARGET :: rivers_data

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ALLOCATE_RIVER_VARS_GRID'

INTEGER ::                                                                     &
   ERROR, error_sum,                                                           &
     ! Error values.
   nx_size, ny_size
     ! Sizes used for allocations.

!end of header
!------------------------------------------------------------------------------

! This routine is only called by the master task. All variables are allocated,
! but they are only allocated at full size when required by the run
! configuration.

! Variables that are always required.
ALLOCATE( grid_riv_pt_number(nx_rivers, ny_rivers), STAT = ERROR )
error_sum = ERROR
ALLOCATE( rivers_dir_int(nx_rivers, ny_rivers), STAT = ERROR )
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rivers_dir(nx_rivers, ny_rivers), STAT = ERROR )
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rivers_lat2d(nx_rivers, ny_rivers), STAT = ERROR )
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rivers_lon2d(nx_rivers, ny_rivers), STAT = ERROR )
error_sum = error_sum + ERROR

! A remapping variable.
ALLOCATE(rivers_data%global_land_index(global_land_pts), STAT = ERROR )
error_sum = error_sum + ERROR

! RFM variables.
IF ( i_river_vn == rivers_rfm ) THEN
  nx_size = nx_rivers
  ny_size = ny_rivers
ELSE
  nx_size = 1
  ny_size = 1
END IF
ALLOCATE(rivers_data%rivers_dra(nx_size,ny_size), STAT = ERROR )
error_sum = error_sum + ERROR

! TRIP variables.
IF ( i_river_vn == rivers_trip ) THEN
  nx_size = nx_rivers
  ny_size = ny_rivers
ELSE
  nx_size = 1
  ny_size = 1
END IF
ALLOCATE(rivers_data%rivers_seq(nx_size,ny_size), STAT = ERROR )
error_sum = error_sum + ERROR

! River coupling variables.
IF ( l_riv_number ) THEN
  nx_size = nx_rivers
  ny_size = ny_rivers
ELSE
  nx_size = 1
  ny_size = 1
END IF
ALLOCATE( rivers_data%rivers_outflow_number(nx_size,ny_size), STAT = ERROR )
error_sum = error_sum + ERROR

IF ( error_sum /= 0 ) THEN
  CALL log_fatal( RoutineName, "Error allocating for rivers arrays." )
END IF

! Initialise to missing data.
grid_riv_pt_number(:,:)       = imdi
rivers_dir_int(:,:)           = imdi
rivers_data%rivers_dir(:,:)   = rmdi
rivers_data%rivers_lat2d(:,:) = rmdi
rivers_data%rivers_lon2d(:,:) = rmdi
rivers_data%global_land_index(:) = imdi
rivers_data%rivers_dra(:,:)   = rmdi
rivers_data%rivers_seq(:,:)   = rmdi
rivers_data%rivers_outflow_number(:,:) = rmdi

! Associate pointers
rivers%rivers_dir   => rivers_data%rivers_dir
rivers%rivers_lat2d => rivers_data%rivers_lat2d
rivers%rivers_lon2d => rivers_data%rivers_lon2d
rivers%global_land_index => rivers_data%global_land_index
rivers%rivers_dra => rivers_data%rivers_dra
rivers%rivers_seq => rivers_data%rivers_seq
rivers%rivers_outflow_number => rivers_data%rivers_outflow_number

RETURN
END SUBROUTINE allocate_river_vars_grid

!##############################################################################

SUBROUTINE process_rivers_data( l_calc_rivers_length, l_use_area,              &
                                rivers, rivers_data )

!------------------------------------------------------------------------------
! Description:
!   Process river ancillary data to produce variables on river points.
!------------------------------------------------------------------------------

USE ancil_info, ONLY: land_pts

USE conversions_mod, ONLY: pi_over_180

USE init_overbank_props_mod, ONLY: allocate_overbank_vars_rp

USE jules_model_environment_mod, ONLY: l_oasis_rivers

USE model_grid_mod, ONLY: global_land_pts, l_coord_latlon,                     &
  projection_x_coord_land, projection_y_coord_land

USE parallel_mod, ONLY: gather_land_field, is_master_task

USE planet_constants_mod, ONLY: planet_radius

USE rivers_regrid_mod, ONLY: calc_map_river_to_land_points

USE jules_rivers_mod, ONLY:                                                    &
  l_riv_overbank, l_trivial_mapping, np_rivers, nx_rivers, ny_rivers,          &
  rivers_dx, rivers_dy, rivers_length, rivers_reglatlon, rivers_regrid,        &
  rivers_x1, rivers_y1,                                                        &
  ! types
  rivers_data_type, rivers_type

USE rivers_utils, ONLY: rivers_earth_area

IMPLICIT NONE

!------------------------------------------------------------------------------
! Scalar arguments with INTENT(IN)
!------------------------------------------------------------------------------
LOGICAL, INTENT(IN) ::                                                         &
  l_calc_rivers_length,                                                        &
    ! Flag indicating if rivers_length will be calculated.
  l_use_area
    ! Switch (used with RFM) to use a drainage area ancillary field to identify
    ! river points.

!------------------------------------------------------------------------------
! Arguments with INTENT(IN OUT)
!------------------------------------------------------------------------------
TYPE(rivers_type), INTENT(IN OUT) :: rivers
TYPE(rivers_data_type), INTENT(IN OUT), TARGET :: rivers_data

!------------------------------------------------------------------------------
! Local parameters.
!------------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'PROCESS_RIVERS_DATA'

!------------------------------------------------------------------------------
! Local scalar variables.
!------------------------------------------------------------------------------
INTEGER ::                                                                     &
  ERROR,                                                                       &
    ! An error value
  error_sum,                                                                   &
    ! Accumulated error values.
  ip, ix, iy, l
    ! Index variables

REAL(KIND=real_jlslsm) ::                                                      &
  half_dlat, half_dlon
    ! Sizes of half a gridbox (degrees)

LOGICAL ::                                                                     &
  l_cyclic_x
    ! Flag indicating a river grid that is cyclic (global) in the x direction.

!------------------------------------------------------------------------------
! Local array variables.
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
  global_proj_x_land(:), global_proj_y_land(:)
    ! Coordinates of land points.

REAL(KIND=real_jlslsm) ::                                                      &
  river_xbound(2), river_ybound(2)
    ! Coordinates of the corners of the river domain - i.e. that part of the
    ! river input grid that will be searched for river points. Units as given
    ! by rivers_xgrid and rivers_ygrid.

!end of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Overview of this subroutine:
! 1. Call process_river_coords to calculate values related to grid coordinates,
!    and the area that will be searched for river points (the river domain).
! 2. Call find_river_points to search the river domain and establish the number
!    of river points.
! 3. Allocate river point arrays.
! 4. Call set_river_point_values to populate the river point arrays using
!    values from the river input grid.
! 5. Call find_downstream_points to identify points downstream of others.
! 6. Calculate river grid mappings.
! Step 3 is done on all tasks, all other steps are only necessary on the
! master task (because rivers are run on that).
!
! Terminology and grids used:
! Land grid:
!   the 2-D grid used by JULES (not rivers) and on which land points sit.
!   If JULES is using a 1-D grid, the model grid is the notional 2D grid
!   across which the points can be scattered.
! River grid:
!   the 2-D river input grid.
! River domain:
!   that part of the river grid that is scanned for river points.
!   The size and location of this depends on the configuration, but ideally it
!   should be just large enough to cover all land points.
! River points:
!    the points within the river domain at which rivers are modelled.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! If l_coord_latlon=T (coordinates are latitude and longitude) this code
! includes provision for reading input files that are presented in N to S order
! before recasting the data in S to N order. Although we could allow the same
! for l_coord_latlon=F, it is simpler to insist that in that case the data are
! presented in strictly increasing order.
!------------------------------------------------------------------------------
IF ( .NOT. l_coord_latlon ) THEN
  DO iy = 2, ny_rivers
    IF ( rivers%rivers_ygrid(iy) <= rivers%rivers_ygrid(iy-1) ) THEN
      CALL log_fatal( RoutineName,                                             &
                      "l_coord_latlon=F requires that river coordinates "   // &
                      "presented in strictly increasing order." )
    END IF
  END DO
END IF

!------------------------------------------------------------------------------
! Gather global land coordinate values from all tasks onto the master task.
!------------------------------------------------------------------------------
ALLOCATE(global_proj_x_land(global_land_pts), STAT = ERROR )
error_sum = ERROR
ALLOCATE(global_proj_y_land(global_land_pts), STAT = ERROR )
error_sum = error_sum + ERROR
IF ( error_sum /= 0 ) THEN
  CALL log_fatal( RoutineName, "Error allocating global lat/lon arrays." )
END IF
! Initialise.
global_proj_x_land(:) = 0.0
global_proj_y_land(:) = 0.0
! Gather values onto master task.
CALL gather_land_field( projection_x_coord_land, global_proj_x_land )
CALL gather_land_field( projection_y_coord_land, global_proj_y_land )

IF ( is_master_task() ) THEN

  !----------------------------------------------------------------------------
  ! Process coordinates and establish the extent of the river domain.
  !----------------------------------------------------------------------------
  CALL process_rivers_coords( global_proj_x_land, global_proj_y_land,          &
             rivers, rivers_dx, rivers_dy, rivers_x1, rivers_y1, l_cyclic_x,   &
             rivers_reglatlon, river_xbound, river_ybound )

  !----------------------------------------------------------------------------
  ! Establish the number and locations of river points.
  !----------------------------------------------------------------------------
  CALL find_river_points( global_proj_x_land,                                  &
                          global_proj_y_land, rivers%rivers_dir,               &
                          rivers%rivers_xgrid, rivers%rivers_ygrid,            &
                          river_xbound, river_ybound, rivers_dir_int,          &
                          np_rivers )

  !----------------------------------------------------------------------------
  ! At this point all values of rivers_dir_int are one of the following:
  !   values in flow_dir_river: these are land values in the river domain
  !   land_off_domain: land but not in the valid river domain
  !   flow_dir_sea:    sea.
  !----------------------------------------------------------------------------

  CALL log_info( RoutineName,                                                  &
                 "River routing points = " // TRIM(to_string(np_rivers)) )
  CALL log_info( RoutineName,                                                  &
                 "Global land pts = " // TRIM(to_string(global_land_pts)) )

  IF ( np_rivers == 0 ) THEN
    CALL log_fatal( RoutineName,                                               &
                    "No river points found. Check inputs." )
  END IF

END IF  ! is_master

!------------------------------------------------------------------------------
! Allocate routing point arrays.
!------------------------------------------------------------------------------
CALL allocate_rivers_vars_rp( np_rivers, rivers, rivers_data )

!------------------------------------------------------------------------------
! Allocate overbank inundation variables.
!------------------------------------------------------------------------------
CALL allocate_overbank_vars_rp( land_pts, np_rivers )

IF ( is_master_task() ) THEN

  !----------------------------------------------------------------------------
  ! Set values at river points.
  !----------------------------------------------------------------------------
  CALL set_river_point_values( l_use_area, rivers_dir_int,                     &
                               rivers, rivers_data, grid_riv_pt_number )

  !----------------------------------------------------------------------------
  ! Identify downstream points.
  !----------------------------------------------------------------------------
  CALL find_downstream_points( l_cyclic_x, rivers%rivers_dir_rp,               &
                               rivers%rivers_index_rp, rivers%rivers_next_rp)

  ! At this point all values of rivers%rivers_next_rp are one of the following:
  !   -10 : an inland drainage point (value = -1*inland_drainage)
  !    -9 : a river mouth (value = -1*river_mouth)
  ! -8:-1 : next point is not a river point - it is outside of the valid river
  !         domain, and possibly off the river grid
  !         (value = -1*rivers%rivers_dir_rp)
  !    >0 : next point is on the grid.

  ! Note that the science code could use rivers%rivers_next_rp values in -8:-1
  ! to sum flow across the edge of the domain, e.g. to check conservation of
  ! mass. That is not currently coded.

  !----------------------------------------------------------------------------
  ! Calculate river length scale and gridbox areas.
  !----------------------------------------------------------------------------
  IF ( l_calc_rivers_length ) THEN
    ! Calculate length scale from grid parameters. Here rivers_dy will be the
    ! gridbox size (degrees), because l_calc_rivers_length can be TRUE only
    ! for lat-lon grids.
    rivers_length = planet_radius * ( ABS( rivers_dy ) * pi_over_180 )
  END IF

  IF ( rivers_reglatlon ) THEN
    half_dlat = 0.5 * rivers_dy
    half_dlon = 0.5 * rivers_dx
    DO ip = 1,np_rivers
      rivers%rivers_boxareas_rp(ip) =                                          &
        ABS( rivers_earth_area( rivers%rivers_lat_rp(ip) - half_dlat,          &
                                rivers%rivers_lat_rp(ip) + half_dlat,          &
                                rivers%rivers_lon_rp(ip) - half_dlon,          &
                                rivers%rivers_lon_rp(ip) + half_dlon ) )
    END DO
  ELSE
    ! Note that this assumes that grids that are not using lat-lon coordinates
    ! are equal area grids - i.e. all gridboxes have the same area.
    rivers%rivers_boxareas_rp(:) = rivers_length * rivers_length
  END IF

  IF ( .NOT. l_oasis_rivers ) THEN
    !--------------------------------------------------------------------------
    ! Calculate mapping ahead of main run.
    !--------------------------------------------------------------------------
    IF ( .NOT. rivers_regrid ) THEN
      ! Find land points that are coincident with river points.
      ! For large grids, this can be slow, hence we issue messages to keep the
      ! user informed.
      CALL log_info( RoutineName,                                              &
                     "Searching for matching river and land points." )
      CALL calc_map_river_to_land_points( global_land_pts,                     &
                                 global_proj_x_land, global_proj_y_land,       &
                                 rivers%rivers_x_coord_rp,                     &
                                 rivers%rivers_y_coord_rp,                     &
                                 rivers%map_river_to_land_points )
      CALL log_info( RoutineName, "Search complete." )

      ! Note that if the mapping between land and river points is trivial
      ! (indicated by l_trivial_mapping=T) we could deallocate
      ! map_river_to_land_points at this point - but that has not beed coded.

    END IF  !  rivers_regrid

    !--------------------------------------------------------------------------
    ! Calculate diagnostics of grid coverage.
    !--------------------------------------------------------------------------
    CALL check_river_coverage( rivers%map_river_to_land_points,                &
                               global_proj_x_land, global_proj_y_land,         &
                               rivers%rivers_x_coord_rp,                       &
                               rivers%rivers_y_coord_rp )

  END IF  !  .NOT. l_oasis_rivers

END IF !  is_master_task

!------------------------------------------------------------------------------
! Dellocate local arrays.
!------------------------------------------------------------------------------
! grid_riv_pt_number is only allocated on the master task, so check.
IF ( ALLOCATED(grid_riv_pt_number) ) THEN
  DEALLOCATE( grid_riv_pt_number )
END IF
DEALLOCATE( global_proj_x_land )
DEALLOCATE( global_proj_y_land )

RETURN
END SUBROUTINE process_rivers_data

!##############################################################################

SUBROUTINE process_rivers_coords( global_proj_x_land, global_proj_y_land,      &
             rivers, rivers_dx, rivers_dy, rivers_x1, rivers_y1, l_cyclic_x,   &
             rivers_reglatlon, river_xbound, river_ybound )

!------------------------------------------------------------------------------
! Description:
!   Process river coordinates, calculating characteristics and checking that
!   grids make sense.
!------------------------------------------------------------------------------

USE jules_rivers_mod, ONLY:                                                    &
  nx_land_grid, nx_rivers, ny_land_grid, ny_rivers, land_dx, land_dy,          &
  rivers_regrid, rivers_type, x1_land_grid, y1_land_grid

USE model_grid_mod, ONLY:                                                      &
  global_land_pts, l_coord_latlon

USE rivers_regrid_mod, ONLY: calc_grid_index

IMPLICIT NONE

!------------------------------------------------------------------------------
! Array arguments with INTENT(IN)
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  global_proj_x_land(global_land_pts), global_proj_y_land(global_land_pts)
    ! Coordinates of land points.

!------------------------------------------------------------------------------
! Type arguments with INTENT(IN OUT)
!------------------------------------------------------------------------------
TYPE(rivers_type), INTENT(IN OUT) :: rivers

!------------------------------------------------------------------------------
! Scalar arguments with INTENT(OUT)
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  rivers_dx,                                                                   &
    ! Size of gridbox of rivers grid in x direction.
    ! If l_coord_latlon=T, the units are degrees of longitude, otherwise units
    ! are unknown.
  rivers_dy,                                                                   &
    ! Size of gridbox of rivers grid in y direction.
    ! If l_coord_latlon=T, the units are degrees of latitude, otherwise units
    ! are unknown.
  rivers_x1,                                                                   &
    ! x coordinate of of westernmost (first) column of gridpoints on the
    ! rivers grid.
  rivers_y1
    ! y coordinate of "southernmost" (first) row of gridpoints on the rivers
    ! grid.

LOGICAL, INTENT(OUT) ::                                                        &
  l_cyclic_x,                                                                  &
    ! Flag indicating a river input grid that is cyclic (global) in the x
    ! direction.
  rivers_reglatlon
    ! Flag indicating if rivers grid is regular in latitude and longitude.

!------------------------------------------------------------------------------
! Array arguments with INTENT(OUT)
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  river_xbound(2), river_ybound(2)
    ! Coordinates of the corners of the river domain - i.e. that part of the
    ! river input grid that will be searched for river points. Units as given
    ! by rivers_xgrid and rivers_ygrid.

!------------------------------------------------------------------------------
! Local parameters.
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm), PARAMETER :: global_lon_threshold = 359.9
  ! Threshold value used when determining if a grid is cyclic (global) in the
  ! longitudinal direction (degrees). This should be close to 360 degrees.

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'PROCESS_RIVERS_COORDS'

!------------------------------------------------------------------------------
! Local scalar variables.
!------------------------------------------------------------------------------
INTEGER :: i, ix, iy, l  ! Index variables

LOGICAL ::                                                                     &
  l_land_180,                                                                  &
    ! Flag indicating the range of land point longitudes.
  l_river_180
    ! Flag indicating the range of river longitudes.

!------------------------------------------------------------------------------
! Local array variables.
!------------------------------------------------------------------------------
INTEGER ::                                                                     &
  land_index_tmp(global_land_pts)
    ! Local version of global land index.

!end of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Report on the range of coordinate values for global land points.
!------------------------------------------------------------------------------
CALL log_info( RoutineName,                                                    &
               'Land points lie in x='                                      // &
               TRIM(to_string( MINVAL(global_proj_x_land) )) // ' to '      // &
               TRIM(to_string( MAXVAL(global_proj_x_land) )) // ' y='       // &
               TRIM(to_string( MINVAL(global_proj_y_land) )) // ' to '      // &
               TRIM(to_string( MAXVAL(global_proj_y_land) )) )

!------------------------------------------------------------------------------
! Get river gridbox size.
!------------------------------------------------------------------------------
rivers_dy = ABS( rivers%rivers_ygrid(2) - rivers%rivers_ygrid(1) )
rivers_dx = ABS( rivers%rivers_xgrid(2) - rivers%rivers_xgrid(1) )
CALL log_info( RoutineName,                                                    &
               'River gridbox size='                                        // &
               TRIM(to_string(rivers_dx)) // ' x '                          // &
               TRIM(to_string(rivers_dy)) )

!------------------------------------------------------------------------------
! Check if the river grid is regular.
!------------------------------------------------------------------------------
! We assume that we always need a regular river grid:
! * if we are regridding, we need regular lat-lon
! * if not regridding (i.e. 1:1 mapping) it is easier if we can use a single
!   length scale to characterise distances and areas.
!
! For each coordinate value we compare with the previous value and flag if
! the separation does not match the expected uniform grid -
! if the difference is a non-negligible fraction of a gridbox.
! Essentially we are checking that this location appears to be on or
! "sufficiently close" to the grid (where "sufficiently close" allows for
! differences from the finite precision of the calculation).
! Note that comparing adjacent values in this way, rather than comparing
! each with an origin, seems to suffer less from issues related to
! precision.

! Check spacings in x are regular.
! We can start at i=3 because 1 and 2 were used to get rivers_dx.
DO i = 3, nx_rivers
  IF ( ABS( ABS( rivers%rivers_xgrid(i) -  rivers%rivers_xgrid(i-1) )          &
            - rivers_dx ) / rivers_dx  > frac_toler ) THEN
    CALL log_fatal( RoutineName, "River input grid must be regular (x)." )
  END IF
END DO

! Check spacings in y are regular.
DO i = 3, ny_rivers
  IF ( ABS( ABS( rivers%rivers_ygrid(i) -  rivers%rivers_ygrid(i-1) )          &
            - rivers_dy ) / rivers_dy  > frac_toler ) THEN
    CALL log_fatal( RoutineName, "River input grid must be regular (y)." )
  END IF
END DO

IF ( l_coord_latlon ) THEN
  ! This is a regular lat-lon river grid.
  rivers_reglatlon = .TRUE.
  CALL log_info( RoutineName, "River grid is a regular lat-lon grid." )
ELSE
  rivers_reglatlon = .FALSE.
END IF

IF ( .NOT. rivers_reglatlon .AND. rivers_regrid ) THEN
  ! Regridding requires a regular lat-lon grid.
  CALL log_fatal( RoutineName,                                                 &
                  "River grid is not a regular lat-lon grid but regridding" // &
                  "is requested. This is not supported." )
ELSE IF ( .NOT. rivers_regrid ) THEN
  ! Remapping requires grids of the same resolution. Here we test that the
  ! difference does not exceed a given fraction of the land resolution.
  IF ( ABS(rivers_dx - land_dx) > frac_toler * land_dx .OR.                    &
       ABS(rivers_dy - land_dy) > frac_toler * land_dy ) THEN
    CALL log_fatal( RoutineName,                                               &
                    "Land and river grids are not of the same resolution " //  &
                    "but regridding is not requested. This is not supported." )
  END IF
END IF

IF ( rivers_reglatlon ) THEN
  !----------------------------------------------------------------------------
  ! We have a regular lat-lon river grid, which requires that the land grid
  ! must also be regular lat-lon (for regridding or mapping). Check that land
  ! points fall on the declared regular grid and return the mapping index.
  !----------------------------------------------------------------------------
  CALL calc_grid_index( global_land_pts, nx_land_grid, ny_land_grid,           &
                        land_dx, land_dy, frac_toler, x1_land_grid,            &
                        y1_land_grid, global_proj_x_land, global_proj_y_land,  &
                        land_index_tmp )
  IF ( rivers_regrid ) THEN
    ! Save the land index.
    rivers%global_land_index(:) = land_index_tmp(:)
  END IF
END IF  !  rivers_reglatlon

!------------------------------------------------------------------------------
! Check that longitudes for land and river points are in the same range.
! This is only required if coordinates are latitude and longitude, because of
! the different ranges that are in common use. For other coordinate systems
! it is assumed that land and river grids use the same general range.
!------------------------------------------------------------------------------
IF ( l_coord_latlon ) THEN
  ! Establish if longitude values of land and river grids are in the same
  ! range. In particular, we want both to either use values in the range
  ! [-180,180] or [0,360]. We assume that the presence of any value < 0
  ! means [-180,180] is used, otherwise it is [0,360].
  IF ( MINVAL( global_proj_x_land ) < 0.0 ) THEN
    l_land_180 = .TRUE.
  ELSE
    l_land_180 = .FALSE.
  END IF
  IF ( MINVAL( rivers%rivers_xgrid ) < 0.0 ) THEN
    l_river_180 = .TRUE.
  ELSE
    l_river_180 = .FALSE.
  END IF
  ! If ranges are different, change river values.
  IF ( l_land_180 .NEQV. l_river_180 ) THEN
    IF ( l_land_180 ) THEN
      ! Land uses [-180,180]. Change river to [-180,180].
      DO ix = 1,nx_rivers
        IF ( rivers%rivers_xgrid(ix) > 180.0 ) THEN
          rivers%rivers_xgrid(ix) = rivers%rivers_xgrid(ix) - 360.0
        END IF
      END DO
    ELSE
      ! l_river_180=T.
      ! Land uses [0,360]. Change river to [0,360].
      DO ix = 1,nx_rivers
        IF ( rivers%rivers_xgrid(ix) < 0.0 ) THEN
          rivers%rivers_xgrid(ix) = rivers%rivers_xgrid(ix) + 360.0
        END IF
      END DO
    END IF
  END IF  !  l_land_180 .NEQV. l_river_180
END IF  !  l_coord_latlon

!------------------------------------------------------------------------------
! Now that we have finished manipulating rivers_xgrid, we can save minima.
! Note that rivers_xgrid might no longer be strictly increasing - there could
! be a discontinuity at 180degE, which is catered for in the code. (It is
! assumed that input coordinates were strictly increasing.)
!------------------------------------------------------------------------------
rivers_x1 = MINVAL( rivers%rivers_xgrid )
rivers_y1 = MINVAL( rivers%rivers_ygrid )

CALL log_info( RoutineName,"River routing input grid: x " //                   &
               TRIM(to_string(MINVAL(rivers%rivers_xgrid)))  //  ' to ' //     &
               TRIM(to_string(MAXVAL(rivers%rivers_xgrid))) //   ' y '  //     &
               TRIM(to_string(MINVAL(rivers%rivers_ygrid))) //   ' to ' //     &
               TRIM(to_string(MAXVAL(rivers%rivers_ygrid))) )

!------------------------------------------------------------------------------
! Decide if the river grid is cyclic (global) in the x direction. We can
! only do this reliably for a regular lat-lon grid - otherwise we have to
! assume that it is not cyclic.
!------------------------------------------------------------------------------
l_cyclic_x = .FALSE.
IF ( rivers_reglatlon ) THEN
  IF ( MAXVAL(rivers%rivers_xgrid) - rivers_x1 + rivers_dx                     &
       > global_lon_threshold ) THEN
    l_cyclic_x = .TRUE.
    CALL log_info( RoutineName, "River grid is cyclic in longitude." )
  END IF
END IF

!------------------------------------------------------------------------------
! Identify the coordinate values that delimit the river domain - i.e. the part
! of the river input grid that is to be searched for river points.
!
! In general the most satisfactory arrangement is if each land gridbox is
! completely covered by river points, and vice versa. In that case all runoff
! (from the land boxes) will be captured by the rivers, and all river boxes
! have well-defined runoff inputs. However:
! (1) Historically JULES could run with river boxes that did not include land
!     (e.g. a large river domain in which is a small land grid is embedded) and
!     we need to be able to reproduce those configurations. The simulated
!     river flow in these river boxes will not be realistic but this
!     configuration could still be considered useful if, e.g., only river
!     points within the land grid are considered and the river network is not
!     missing upstream areas that contribute a significant amount of runoff.
! (2) In general it is non-trivial to reliably establish the overlap between
!     grids which can be of different resolutions.
!------------------------------------------------------------------------------
!
! Use values from the namelist. The coordinates used to outline the valid
! domain are calculated as in older versions - essentially assuming that input
! parameters (e.g. y1_land_grid, ny_land_grid) are reasonable.
! If these inputs are not specified correctly, the model might still run but the
! results can be dubious - e.g. river points can be set up outside the area of
! the land grid but these will not have any runoff inputs, meaning the modelled
! riverflows are compromised. In other instances an error will be raised at a
! later stage.
!
! In the older code, river gridpoints are simply compared with these
! notional land extremes - and we do that here.
! If we assume that x1_land_grid and y1_land_grid have been set to indicate the
! gridpoint in the SW corner of the land grid, and nx_land_grid and ny_land_grid
! are the smallest values that include all land points, then by using
! nx_land_grid (rather than nx_land_grid-1) this calculation effectively allows
! an extra column of points at the E edge (and similarly in the y direction),
! i.e. the domain would be larger than required.
river_xbound(1) = x1_land_grid
river_xbound(2) = x1_land_grid + nx_land_grid * land_dx
river_ybound(1) = y1_land_grid
river_ybound(2) = y1_land_grid + ny_land_grid * land_dy

CALL log_info( RoutineName,                                                    &
               "River points will be searched for in x "                    // &
               TRIM(to_string(river_xbound(1))) // ' to '                   // &
               TRIM(to_string(river_xbound(2))) // ' y '                    // &
               TRIM(to_string(river_ybound(1))) // ' to '                   // &
               TRIM(to_string(river_ybound(2))) )

RETURN
END SUBROUTINE process_rivers_coords

!##############################################################################

SUBROUTINE find_river_points( global_proj_x_land,                              &
                              global_proj_y_land, rivers_dir, rivers_xgrid,    &
                              rivers_ygrid, river_xbound, river_ybound,        &
                              rivers_dir_int, np_rivers )

!------------------------------------------------------------------------------
! Description:
!   Identify points that are rivers.
!------------------------------------------------------------------------------

USE jules_rivers_mod, ONLY:                                                    &
  nx_rivers, ny_rivers, rivers_regrid

USE model_grid_mod, ONLY: global_land_pts

IMPLICIT NONE

!------------------------------------------------------------------------------
! Array arguments with INTENT(IN)
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  global_proj_x_land(global_land_pts), global_proj_y_land(global_land_pts),    &
    ! Coordinates of land points.
  rivers_dir(nx_rivers,ny_rivers),                                             &
    ! River routing direction index.
  river_xbound(2), river_ybound(2),                                            &
    ! Coordinates of the corners of the river domain - i.e. that part of the
    ! river input grid that will be searched for river points. Units as given
    ! by rivers_xgrid and rivers_ygrid.
  rivers_xgrid(nx_rivers), rivers_ygrid(ny_rivers)
    ! Coordinates for rivers grid

!------------------------------------------------------------------------------
! Array arguments with INTENT(IN OUT)
!------------------------------------------------------------------------------
INTEGER, INTENT(IN OUT) ::                                                     &
  rivers_dir_int(nx_rivers, ny_rivers)
    ! River routing direction index (integer version).

!------------------------------------------------------------------------------
! Scalar arguments with INTENT(OUT)
!------------------------------------------------------------------------------
INTEGER, INTENT(OUT) ::                                                        &
  np_rivers
    ! Number of river points.

!------------------------------------------------------------------------------
! Local parameters.
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm), PARAMETER :: coord_toler = 1.0e-6
    ! Tolerance allowed when comparing coordinates.
    ! If l_coord_latlon=T, the units are degrees of latitude and longitude,
    ! otherwise units are unknown.

INTEGER ::                                                                     &
  ix, iy, l
    ! Loop counters.

REAL(KIND=real_jlslsm) ::                                                      &
  x_coord, y_coord
    ! Coordinate values.

LOGICAL ::                                                                     &
  l_match
    ! Used when identifying matching points.

!end of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Scan the grid looking for river points.
!------------------------------------------------------------------------------
np_rivers = 0

DO iy = 1,ny_rivers
  DO ix = 1,nx_rivers

    !--------------------------------------------------------------------------
    ! Create an INTEGER version of rivers_dir. The REAL version is not
    ! used after this.
    !--------------------------------------------------------------------------
    rivers_dir_int(ix,iy) = NINT( rivers_dir(ix,iy) )

    !-------------------------------------------------------------------------
    ! Consider points on the valid domain only (i.e. the part of the river
    ! grid that is to be searched for river points).
    !-------------------------------------------------------------------------
    IF ( ( rivers_xgrid(ix) < river_xbound(1) ) .OR.                           &
         ( rivers_ygrid(iy) < river_ybound(1) ) .OR.                           &
         ( rivers_xgrid(ix) > river_xbound(2) ) .OR.                           &
         ( rivers_ygrid(iy) > river_ybound(2) ) ) THEN

      !-----------------------------------------------------------------------
      ! This point is not in the valid river domain.
      ! If this is a land point (as indicated by the flow direction) set
      ! the direction to land_off_domain, or if this is sea (all other
      ! values) set to flow_dir_sea.
      !-----------------------------------------------------------------------
      IF ( ANY( rivers_dir_int(ix,iy) == flow_dir_river(:) ) ) THEN
        rivers_dir_int(ix,iy) = land_off_domain
      ELSE
        rivers_dir_int(ix,iy) = flow_dir_sea
      END IF

    ELSE IF ( ANY( rivers_dir_int(ix,iy) == flow_dir_river(:) ) ) THEN

      !----------------------------------------------------------------------
      ! This point is in the valid river domain and is land (according to
      ! rivers_dir_int). It is a river point.
      !----------------------------------------------------------------------
      np_rivers = np_rivers + 1

    ELSE

      !----------------------------------------------------------------------
      ! This point is in the valid river domain but is not land (according to
      ! flow_dir_int).
      ! Set flow direction to indicate sea. This is not currently necessary
      ! but doing so ensures that the same value is used for sea across
      ! the grid, in case any later code changes rely on that.
      !----------------------------------------------------------------------
      rivers_dir_int(ix,iy) = flow_dir_sea

    END IF

  END DO  !  ix
END DO  !  iy

RETURN
END SUBROUTINE find_river_points

!#############################################################################

SUBROUTINE allocate_rivers_vars_rp( np_rivers, rivers, rivers_data )

!------------------------------------------------------------------------------
! Description:
!   Allocate river point variables, initialise, and associate pointers.
!------------------------------------------------------------------------------

USE jules_model_environment_mod, ONLY: l_oasis_rivers

USE jules_rivers_mod, ONLY: i_river_vn, rivers_data_type, rivers_rfm,          &
                            rivers_trip, rivers_type

USE parallel_mod, ONLY: is_master_task

IMPLICIT NONE

!------------------------------------------------------------------------------
! Scalar arguments with INTENT(IN)
!------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  np_rivers
    ! Number of river points.

!------------------------------------------------------------------------------
! Arguments with INTENT(IN OUT)
!------------------------------------------------------------------------------
TYPE(rivers_type), INTENT(IN OUT) :: rivers
TYPE(rivers_data_type), INTENT(IN OUT), TARGET :: rivers_data

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ALLOCATE_RIVERS_VARS_RP'

!------------------------------------------------------------------------------
! Local scalar variables.
!------------------------------------------------------------------------------
INTEGER ::                                                                     &
  ERROR, error_sum,                                                            &
    ! Error flags.
  nland_tmp,                                                                   &
    ! Number of land points to allocate for.
  np_rivers_tmp
    ! Number of river points to allocate for.

!end of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! As river routing is performed on the master task, we only allocate arrays at
! full size on that task, and only if the run configuration requires it.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Allocate variables that are required by all configurations.
!------------------------------------------------------------------------------
IF ( is_master_task() ) THEN
  ! Full size.
  np_rivers_tmp = np_rivers
ELSE
  ! Minimum size.
  np_rivers_tmp = 1
END IF

! Ancillaries and related fields.
ALLOCATE(rivers_data%rivers_boxareas_rp(np_rivers_tmp),  STAT = ERROR)
error_sum = ERROR
ALLOCATE(rivers_data%rivers_dir_rp(np_rivers_tmp),       STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rivers_index_rp(np_rivers_tmp),     STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rivers_lat_rp(np_rivers_tmp),       STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rivers_lon_rp(np_rivers_tmp),       STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rivers_next_rp(np_rivers_tmp),      STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rivers_x_coord_rp(np_rivers_tmp),   STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rivers_y_coord_rp(np_rivers_tmp),   STAT = ERROR)
error_sum = error_sum + ERROR
! Fluxes.
ALLOCATE(rivers_data%rflow_rp(np_rivers_tmp),            STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%riverout_rp(np_rivers),   STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rrun_rp(np_rivers_tmp),             STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rrun_sub_surf_rp(np_rivers_tmp),    STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rrun_surf_rp(np_rivers_tmp),        STAT = ERROR)
error_sum = error_sum + ERROR

!------------------------------------------------------------------------------
! Allocate variable specific to remapping.
!------------------------------------------------------------------------------
ALLOCATE(rivers_data%map_river_to_land_points(np_rivers_tmp), STAT = ERROR)
error_sum = error_sum + ERROR

!------------------------------------------------------------------------------
! Allocate RFM variables
!------------------------------------------------------------------------------
IF ( i_river_vn == rivers_rfm .AND. is_master_task() ) THEN
  ! Full size.
  np_rivers_tmp = np_rivers
ELSE
  ! Minimum size.
  np_rivers_tmp = 1
END IF

! Ancillaries and related fields.
! Note: rfm_flowobs1_rp is not used by standalone JULES - but is left here
!       for now.
ALLOCATE(rivers_data%rfm_flowobs1_rp(np_rivers_tmp),  STAT = ERROR )
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rfm_iarea_rp(np_rivers_tmp),     STAT = ERROR )
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rfm_land_rp(np_rivers_tmp),      STAT = ERROR )
error_sum = error_sum + ERROR
! Fluxes and stores.
ALLOCATE(rivers_data%rfm_baseflow_rp(np_rivers_tmp),  STAT = ERROR )
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rfm_bflowin_rp(np_rivers_tmp),   STAT = ERROR )
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rfm_flowin_rp(np_rivers_tmp),    STAT = ERROR )
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rfm_rivflow_rp(np_rivers_tmp),   STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rfm_substore_rp(np_rivers_tmp),  STAT = ERROR )
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rfm_surfstore_rp(np_rivers_tmp), STAT = ERROR )
error_sum = error_sum + ERROR

!------------------------------------------------------------------------------
! Allocate TRIP variables
!------------------------------------------------------------------------------
IF ( i_river_vn == rivers_trip .AND. is_master_task() ) THEN
  ! Full size.
  np_rivers_tmp = np_rivers
ELSE
  ! Minimum size.
  np_rivers_tmp = 1
END IF
! Ancillaries and related fields.
ALLOCATE(rivers_data%rivers_seq_rp(np_rivers_tmp), STAT = ERROR)
error_sum = error_sum + ERROR
! Stores.
ALLOCATE(rivers_data%rivers_sto_rp(np_rivers_tmp), STAT = ERROR)
error_sum = error_sum + ERROR

!------------------------------------------------------------------------------
! Allocate coupling variables.
!------------------------------------------------------------------------------
IF ( l_oasis_rivers .AND. is_master_task() ) THEN
  ! Full size.
  np_rivers_tmp = np_rivers
ELSE
  ! Minimum size.
  np_rivers_tmp = 1
END IF
ALLOCATE(rivers_data%sub_surf_roff_rp(np_rivers_tmp), STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%surf_roff_rp(np_rivers_tmp), STAT = ERROR)
error_sum = error_sum + ERROR
ALLOCATE(rivers_data%rivers_outflow_number_rp(np_rivers_tmp), STAT = ERROR)
error_sum = error_sum + ERROR

IF ( error_sum /= 0 ) THEN
  CALL log_fatal( RoutineName, "Error allocating for routing point arrays." )
END IF

!------------------------------------------------------------------------------
! Initialise array values to missing values, except fluxes can be zero to
! preserve existing results.
!------------------------------------------------------------------------------
! Ancillaries etc.
rivers_data%rivers_boxareas_rp(:) = rmdi
rivers_data%rivers_dir_rp(:)     = imdi
rivers_data%rivers_index_rp(:)   = imdi
rivers_data%rivers_lat_rp(:)     = rmdi
rivers_data%rivers_lon_rp(:)     = rmdi
rivers_data%rivers_next_rp(:)    = imdi
rivers_data%rivers_x_coord_rp(:) = rmdi
rivers_data%rivers_y_coord_rp(:) = rmdi
! Fluxes.
rivers_data%rflow_rp(:)          = 0.0
IF ( l_oasis_rivers ) THEN
  ! For now preserve coupling field
  rivers_data%riverout_rp(:)     = 0.0
ELSE
  rivers_data%riverout_rp(:)     = rmdi
END IF
rivers_data%rrun_rp(:)           = 0.0
rivers_data%rrun_sub_surf_rp(:)  = 0.0
rivers_data%rrun_surf_rp(:)      = 0.0

! Initialise regridding variable.
rivers_data%map_river_to_land_points(:) = imdi

! Initialise RFM array values
rivers_data%rfm_flowobs1_rp(:)  = rmdi
rivers_data%rfm_iarea_rp(:)     = imdi
rivers_data%rfm_land_rp(:)      = imdi
rivers_data%rfm_baseflow_rp(:)  = 0.0
! rfm_bflowin_rp, rfm_flowin_rp, rfm_substore_rp and rfm_surfstore_rp should
! be initiailsed to zero to preserve existing results and until better
! initialisation is possible.
rivers_data%rfm_bflowin_rp(:)   = 0.0
rivers_data%rfm_flowin_rp(:)    = 0.0
rivers_data%rfm_rivflow_rp(:)   = 0.0
rivers_data%rfm_substore_rp(:)  = 0.0
rivers_data%rfm_surfstore_rp(:) = 0.0

! Initialise TRIP array values
rivers_data%rivers_seq_rp(:) = imdi
! rivers_sto_rp should be initialised to zero to preserve existing results and
! until better initialisation is possible.
rivers_data%rivers_sto_rp(:) = 0.0

! Initialise coupling variables.
rivers_data%rivers_outflow_number_rp(:) = imdi
rivers_data%sub_surf_roff_rp(:) = 0.0
rivers_data%surf_roff_rp(:)     = 0.0

!------------------------------------------------------------------------------
! Associate pointers
!------------------------------------------------------------------------------
! Ancillaries etc.
rivers%rivers_boxareas_rp => rivers_data%rivers_boxareas_rp
rivers%rivers_dir_rp => rivers_data%rivers_dir_rp
rivers%rivers_index_rp => rivers_data%rivers_index_rp
rivers%rivers_lat_rp => rivers_data%rivers_lat_rp
rivers%rivers_lon_rp => rivers_data%rivers_lon_rp
rivers%rivers_next_rp => rivers_data%rivers_next_rp
rivers%rivers_x_coord_rp => rivers_data%rivers_x_coord_rp
rivers%rivers_y_coord_rp => rivers_data%rivers_y_coord_rp
! Fluxes.
rivers%rflow_rp => rivers_data%rflow_rp
rivers%riverout_rp => rivers_data%riverout_rp
rivers%rrun_rp => rivers_data%rrun_rp
rivers%rrun_sub_surf_rp => rivers_data%rrun_sub_surf_rp
rivers%rrun_surf_rp => rivers_data%rrun_surf_rp

! Associate pointers for regridding variable.
rivers%map_river_to_land_points => rivers_data%map_river_to_land_points

! Associate pointers for RFM variables.
rivers%rfm_flowobs1_rp => rivers_data%rfm_flowobs1_rp
rivers%rfm_iarea_rp => rivers_data%rfm_iarea_rp
rivers%rfm_land_rp => rivers_data%rfm_land_rp
rivers%rfm_baseflow_rp => rivers_data%rfm_baseflow_rp
rivers%rfm_bflowin_rp => rivers_data%rfm_bflowin_rp
rivers%rfm_flowin_rp => rivers_data%rfm_flowin_rp
rivers%rfm_rivflow_rp => rivers_data%rfm_rivflow_rp
rivers%rfm_substore_rp => rivers_data%rfm_substore_rp
rivers%rfm_surfstore_rp => rivers_data%rfm_surfstore_rp

! Associate pointers for TRIP variables.
rivers%rivers_seq_rp => rivers_data%rivers_seq_rp
rivers%rivers_sto_rp => rivers_data%rivers_sto_rp

! Associate pointers for coupling variables.
rivers%rivers_outflow_number_rp => rivers_data%rivers_outflow_number_rp
rivers%sub_surf_roff_rp => rivers_data%sub_surf_roff_rp
rivers%surf_roff_rp => rivers_data%surf_roff_rp

RETURN
END SUBROUTINE allocate_rivers_vars_rp

!##############################################################################

SUBROUTINE set_river_point_values( l_use_area, rivers_dir_int,                 &
                                   rivers, rivers_data, grid_riv_pt_number )

!------------------------------------------------------------------------------
! Description:
!   Set values at river points - generally by extracting from 2-D ancillary
!   fields.
!------------------------------------------------------------------------------

USE jules_model_environment_mod, ONLY: l_oasis_rivers

USE jules_rivers_mod, ONLY:                                                    &
  a_thresh, i_river_vn, l_riv_overbank, np_rivers, nseqmax, nx_rivers,         &
  ny_rivers, rfm_land, rfm_river, rivers_dx, rivers_dy, rivers_rfm,            &
  rivers_trip, rivers_x1,                                                      &
  ! types
  rivers_data_type, rivers_type

USE model_grid_mod, ONLY: l_coord_latlon

USE overbank_inundation_mod, ONLY:                                             &
  logn_mean, logn_stdev, logn_mean_rp, logn_stdev_rp, overbank_hypsometric,    &
  overbank_model

IMPLICIT NONE

!------------------------------------------------------------------------------
! Scalar arguments with INTENT(IN)
!------------------------------------------------------------------------------
LOGICAL, INTENT(IN) ::                                                         &
  l_use_area
    ! Switch (used with RFM) to use a drainage area ancillary field to identify
    ! river points.

!------------------------------------------------------------------------------
! Array arguments with INTENT(IN)
!------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  rivers_dir_int(nx_rivers, ny_rivers)
    ! Integer version of rivers_dir.

!------------------------------------------------------------------------------
! Arguments with INTENT(IN OUT)
!------------------------------------------------------------------------------
TYPE(rivers_type), INTENT(IN OUT) :: rivers
TYPE(rivers_data_type), INTENT(IN OUT), TARGET :: rivers_data

!------------------------------------------------------------------------------
! Array arguments with INTENT(OUT)
!------------------------------------------------------------------------------
INTEGER, INTENT(OUT) ::                                                        &
  grid_riv_pt_number(nx_rivers, ny_rivers)
    ! Map full river grid to river points. For each location on the grid, this
    ! is the river point number.

!------------------------------------------------------------------------------
! Local parameters.
!------------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SET_RIVER_POINT_VALUES'

!------------------------------------------------------------------------------
! Local scalar variables.
!------------------------------------------------------------------------------
INTEGER ::                                                                     &
  ERROR, error_sum,                                                            &
    ! Error variables.
  i1, i2,                                                                      &
    ! Loop bounds.
  ip, irx, iry, ix, iy,                                                        &
    ! Counters and index variables.
  step
    ! Loop increment.

!end of header
!------------------------------------------------------------------------------

! Initialise mapping to show no points on grid.
grid_riv_pt_number(:,:) = 0

!----------------------------------------------------------------------------
! Ensure that the river grid is scanned from "S to N" (meaning the coordinate
! value increases from y=1 to y=ny_rivers; this coordinate might not be
! latitude) irrespective of how values are stored in the ancillary files.
!----------------------------------------------------------------------------
IF ( rivers%rivers_ygrid(1) > MINVAL(rivers%rivers_ygrid) ) THEN
  i1   = ny_rivers
  i2   = 1
  step = -1
ELSE
  i1   = 1
  i2   = ny_rivers
  step = 1
END IF

! Initialise max sequence and point number to zero.
nseqmax = 0
ip      = 0

!----------------------------------------------------------------------------
! Here we are scanning the river input grid, with indices ix and iy.
! If point (ix,iy) is a river point, we populate the river grid at location
! (irx,iry).
! irx is calculated so as to ensure that the river grid is formed with a
! strictly increasing x coordinate - taking account of any adjustment that
! was made to match the range of longitude to that of the land grid.
! iry is calculated so as to ensure that the river grid is formed with a
! strictly increasing y coordinate.
!----------------------------------------------------------------------------
DO ix = 1,nx_rivers

  ! Calculate x position on river grid.
  ! Note that if l_coord_latlon=F or we have not adjusted longitude values,
  ! (i.e. rivers_x1 = rivers%rivers_xgrid(1)) this will result in irx=1,2,...
  ! This ensures that the river grid (the river domain) has a strictly
  ! increasing x coordinate.
  irx = NINT( (rivers%rivers_xgrid(ix) - rivers_x1) / rivers_dx ) + 1

  ! Initialise y position on river grid.
  iry = 0

  ! This loop over y will scan up or down columns, whichever way is required
  ! to ensure the coordinate value is increasing bottom to top on the river
  ! domain.
  DO iy = i1,i2,step

    ! We fill the domain bottom to top.
    iry = iry + 1

    ! Populate 2-D latitude and longitude coordinates.
    ! Note that if l_coord_latlon=F, the lat-lon coordinates were read from
    ! a file and there is nothing to do here.
    IF ( l_coord_latlon ) THEN
      rivers%rivers_lat2d(irx,iry) = rivers%rivers_ygrid(iy)
      rivers%rivers_lon2d(irx,iry) = rivers%rivers_xgrid(ix)
    END IF

    !------------------------------------------------------------------------
    ! Set values at river points.
    !------------------------------------------------------------------------
    IF ( ANY( rivers_dir_int(ix,iy) == flow_dir_river(:) ) ) THEN

      ip = ip + 1


      !----------------------------------------------------------------------
      ! Define mapping between river vector and routing grid.
      !----------------------------------------------------------------------
      grid_riv_pt_number(irx,iry) = ip

      !----------------------------------------------------------------------
      ! Calculate index of position on river input grid.
      ! Note that this index scans row-wise from bottom left, bottom to top.
      ! The loops that control the order in which points are found do not
      ! necessarily scan the same way - which is confusing but acceptable.
      !----------------------------------------------------------------------
      rivers%rivers_index_rp(ip) = (iry-1) * nx_rivers + irx

      !----------------------------------------------------------------------
      ! Set river vector information based on 2D input ancillaries.
      ! First, coordinates and flow direction.
      !----------------------------------------------------------------------
      rivers%rivers_lon_rp(ip)     = rivers%rivers_lon2d(irx,iry)
      rivers%rivers_lat_rp(ip)     = rivers%rivers_lat2d(irx,iry)
      rivers%rivers_x_coord_rp(ip) = rivers%rivers_xgrid(ix)
      rivers%rivers_y_coord_rp(ip) = rivers%rivers_ygrid(iy)
      rivers%rivers_dir_rp(ip)     = rivers_dir_int(ix,iy)

      !----------------------------------------------------------------------
      ! Set values that depend on the chosen river model.
      !----------------------------------------------------------------------
      IF ( i_river_vn == rivers_rfm ) THEN
        !--------------------------------------------------------------------
        ! RFM variables.
        !--------------------------------------------------------------------
        IF ( l_use_area ) THEN
          ! Set drainage area from ancillary.
          rivers%rfm_iarea_rp(ip) = NINT( rivers%rivers_dra(ix,iy) )
        ELSE
          ! Ancillary not present. Set area above threshold, so that all
          ! points will be identified as rivers.
          rivers%rfm_iarea_rp(ip) = a_thresh + 1
        END IF
        ! Set points to be either land or river.
        IF (rivers%rfm_iarea_rp(ip) > a_thresh) THEN
          rivers%rfm_land_rp(ip) = rfm_river  ! river
        ELSE
          rivers%rfm_land_rp(ip) = rfm_land   ! land
        END IF

      ELSE IF ( i_river_vn == rivers_trip ) THEN
        !--------------------------------------------------------------------
        ! TRIP variables.
        !--------------------------------------------------------------------
        rivers%rivers_seq_rp(ip) = NINT(rivers%rivers_seq(ix,iy))
        IF ( rivers%rivers_seq_rp(ip) > nseqmax ) THEN
          nseqmax = rivers%rivers_seq_rp(ip)
        END IF

      END IF  !  i_river_vn

      !----------------------------------------------------------------------
      ! Set values of coupling variables.
      !----------------------------------------------------------------------
      IF ( l_oasis_rivers ) THEN
        rivers%rivers_outflow_number_rp(ip) =                                  &
                                    NINT(rivers%rivers_outflow_number(ix,iy))
      END IF

      !----------------------------------------------------------------------
      ! Set overbank inundation ancillary variables, if required.
      !----------------------------------------------------------------------
      IF ( l_riv_overbank .AND. overbank_model == overbank_hypsometric ) THEN
        logn_mean_rp(ip)  = logn_mean(ix,iy)
        logn_stdev_rp(ip) = logn_stdev(ix,iy)
      END IF

    END IF  !  rivers_dir_int

  END DO  !  iy river rows
END DO  !  ix river columns

! Double check that we have found all the river points.
IF ( ip /= np_rivers ) THEN
  CALL log_fatal( RoutineName,                                                 &
                  "Failure of logic; did not find all river points." )
END IF

RETURN
END SUBROUTINE set_river_point_values

!##############################################################################

SUBROUTINE find_downstream_points( l_cyclic_x, rivers_dir_rp, rivers_index_rp, &
                                   rivers_next_rp )


!------------------------------------------------------------------------------
! Description:
!   Identify the next downstream point from flow direction.
!------------------------------------------------------------------------------

USE jules_rivers_mod, ONLY:                                                    &
  flow_dir_delta, np_rivers, nx_rivers, ny_rivers

USE rivers_regrid_mod, ONLY: get_xy_pos

IMPLICIT NONE

!------------------------------------------------------------------------------
! Scalar arguments with INTENT(IN)
!------------------------------------------------------------------------------
LOGICAL, INTENT(IN) ::                                                         &
  l_cyclic_x
    ! Flag indicating a river grid that is cyclic (global) in the x direction.

!------------------------------------------------------------------------------
! Array arguments with INTENT(IN)
!------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  rivers_dir_rp(np_rivers),                                                    &
    ! River routing direction index
  rivers_index_rp(np_rivers)

!------------------------------------------------------------------------------
! Array arguments with INTENT(OUT)
!------------------------------------------------------------------------------
INTEGER, INTENT(OUT) ::                                                        &
  rivers_next_rp(np_rivers)

!------------------------------------------------------------------------------
! Local parameters.
!------------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'FIND_DOWNSTREAM_POINTS'

!------------------------------------------------------------------------------
! Local scalar variables.
!------------------------------------------------------------------------------
INTEGER ::                                                                     &
  inext, ip, ix, iy, jnext
    ! Counters and index variables.

!end of header
!------------------------------------------------------------------------------

! Note that there is no need to initialise rivers_next_rp as we set a value at
! every point.

DO ip = 1,np_rivers

  SELECT CASE ( rivers_dir_rp(ip) )

  CASE ( 1:8 )
    ! These are the 8 "physical" flow directions, i.e. there is a flow
    ! direction and a point downstream of this.
    CALL get_xy_pos(rivers_index_rp(ip),nx_rivers,ny_rivers,ix,iy)
    inext = ix + flow_dir_delta( rivers_dir_rp(ip),1 )
    jnext = iy + flow_dir_delta( rivers_dir_rp(ip),2 )

    ! Deal with cyclic boundaries.
    IF ( l_cyclic_x ) THEN
      IF ( inext > nx_rivers ) inext = inext - nx_rivers
      IF ( inext < 1 ) inext = inext + nx_rivers
    END IF

    ! Identify the next point downstream.
    IF ( inext >= 1 .AND. inext <= nx_rivers .AND. jnext >=1                   &
         .AND. jnext <= ny_rivers ) THEN
      ! The next point is on the river grid. Check that it is a river point.
      IF ( grid_riv_pt_number(inext, jnext) > 0 ) THEN
        ! The next point is a river point.
        rivers_next_rp(ip) = grid_riv_pt_number(inext,jnext)
      ELSE
        ! The next point is on the grid but is not a river point.
        IF ( rivers_dir_int(inext,jnext) == flow_dir_sea ) THEN
          ! The next point is sea. In many ancillary files this point would
          ! have been flagged as a river mouth. Convert the current  point
          ! to a river mouth (but there is no need to change
          ! rivers_dir_rp as it is not used again).
          rivers_next_rp(ip) = -river_mouth
        ELSE
          ! The next point is land (but it is not in the valid river
          ! domain). Save the direction in which it lies.
          ! Note that because it it on the grid we know its location (i.e.
          ! we have lat-lon values) but this information will not be
          ! available to the science code (which only gets the coordinates
          ! of river points).
          rivers_next_rp(ip) = -rivers_dir_rp(ip)
        END IF
      END IF  !  grid_riv_pt_number
    ELSE
      ! The next point is not on the grid, but we know which direction it
      ! lies in. Save that information for later use when getting a length
      ! for the flow. Note that flow across the E or W edge of a cyclic
      ! global grid is dealt with above.
      rivers_next_rp(ip) = -rivers_dir_rp(ip)
    END IF  !  inext

  CASE ( river_mouth )

    ! Identify this with a special value < 0.
    rivers_next_rp(ip) = -river_mouth

  CASE ( inland_drainage )
    ! Identify this with a special value < 0.
    rivers_next_rp(ip) = -inland_drainage

  CASE DEFAULT
    CALL log_fatal( RoutineName,                                               &
                    "Unexpected flow direction: "                           // &
                    TRIM(to_string(rivers_dir_rp(ip))) )

  END SELECT  !  rivers_dir_rp

END DO  !  ip

END SUBROUTINE find_downstream_points

!##############################################################################

SUBROUTINE check_river_coverage( map_river_to_land_points, global_proj_x_land, &
                 global_proj_y_land, rivers_x_coord_rp, rivers_y_coord_rp )

!------------------------------------------------------------------------------
! Description:
!   Assesses the extent to which the land and river layers overlap, so as to
!   produce diagnostic messages.
!------------------------------------------------------------------------------

USE jules_rivers_mod, ONLY: land_dx, land_dy, l_trivial_mapping, np_rivers,    &
                            rivers_dx, rivers_dy, rivers_regrid

USE model_grid_mod, ONLY: global_land_pts

IMPLICIT NONE

!------------------------------------------------------------------------------
! Array arguments with INTENT(IN)
!------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  map_river_to_land_points(np_rivers)
    ! List of coincident land point numbers, on river points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  global_proj_x_land(global_land_pts),                                         &
    ! x coordinates of land points.
  global_proj_y_land(global_land_pts),                                         &
    ! y coordinates of land points.
  rivers_x_coord_rp(np_rivers),                                                &
    ! River routing point projection x coordinates.
  rivers_y_coord_rp(np_rivers)
    ! River routing point projection y coordinates.

!------------------------------------------------------------------------------
! Local scalar parameters.
!------------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CHECK_RIVER_COVERAGE'

!------------------------------------------------------------------------------
! Local scalar variables.
!------------------------------------------------------------------------------
INTEGER ::                                                                     &
  num_match,                                                                   &
    ! Number of points for which a matching point is found.
  num_no_match
    ! Number of points for which no matching point is found.

REAL(KIND=real_jlslsm) ::                                                      &
  dx_toler, dy_toler
    ! Tolerances. Units depend on the coordinate system in use.

!------------------------------------------------------------------------------
! Local array variables.
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  land_x_bounds_edge(2),                                                       &
    ! x coordinates of westermost and easternmost edges of land gridboxes.
  land_y_bounds_edge(2),                                                       &
    ! y coordinates of southernmost and northernmost edges of land gridboxes.
  rivers_x_bounds_edge(2),                                                     &
    ! x coordinates of westermost and easternmost edges of river gridboxes.
  rivers_y_bounds_edge(2)
    ! y coordinates of southernmost and northernmost edges of river gridboxes.

!end of header
!------------------------------------------------------------------------------

IF ( rivers_regrid ) THEN

  !----------------------------------------------------------------------------
  ! Regridding can take place between grids of different resolutions and with
  ! grids that do not nest (i.e. a gridbox on one grid can straddle multiple
  ! gridboxes on the other grid). As such it is non-trivial to determine the
  ! extent to which each land gridbox is covered by river boxes, and vice
  ! versa. Instead we test and report on some necessary but not sufficient
  ! conditions.
  !----------------------------------------------------------------------------

  ! Get the coordinates at the outer edges of the extreme land gridboxes.
  land_x_bounds_edge(1) = MINVAL( global_proj_x_land ) - 0.5 * land_dx
  land_x_bounds_edge(2) = MAXVAL( global_proj_x_land ) + 0.5 * land_dx
  land_y_bounds_edge(1) = MINVAL( global_proj_y_land ) - 0.5 * land_dy
  land_y_bounds_edge(2) = MAXVAL( global_proj_y_land ) + 0.5 * land_dy

  CALL log_info( RoutineName,                                                  &
                 "Outer edges of land gridboxes are at x="                 //  &
                  TRIM(to_string(land_x_bounds_edge(1))) // ' to '         //  &
                  TRIM(to_string(land_x_bounds_edge(2))) // ' y='          //  &
                  TRIM(to_string(land_y_bounds_edge(1))) // ' to '         //  &
                  TRIM(to_string(land_y_bounds_edge(2))) )

  ! Get the coordinates at the outer edges of the extreme river gridboxes.
  rivers_x_bounds_edge(1) = MINVAL( rivers_x_coord_rp(:) ) - 0.5 * rivers_dx
  rivers_x_bounds_edge(2) = MAXVAL( rivers_x_coord_rp(:) ) + 0.5 * rivers_dx
  rivers_y_bounds_edge(1) = MINVAL( rivers_y_coord_rp(:) ) - 0.5 * rivers_dy
  rivers_y_bounds_edge(2) = MAXVAL( rivers_y_coord_rp(:) ) + 0.5 * rivers_dy

  CALL log_info( RoutineName,                                                  &
                 "Outer edges of river gridboxes are at x="                 // &
                 TRIM(to_string(rivers_x_bounds_edge(1))) // ' to '         // &
                 TRIM(to_string(rivers_x_bounds_edge(2))) // ' y='          // &
                 TRIM(to_string(rivers_y_bounds_edge(1))) // ' to '         // &
                 TRIM(to_string(rivers_y_bounds_edge(2))) )

  ! Look for river gridboxes that have some area outside the land area, by
  ! comparing the coordinates of the edges of the extreme gridboxes, allowing
  ! a small tolerance based on river gridbox size.
  dx_toler = frac_toler * rivers_dx
  dy_toler = frac_toler * rivers_dy
  IF ( rivers_x_bounds_edge(1) < land_x_bounds_edge(1) - dx_toler  .OR.        &
       rivers_x_bounds_edge(2) > land_x_bounds_edge(2) + dx_toler  .OR.        &
       rivers_y_bounds_edge(1) < land_y_bounds_edge(1) - dy_toler  .OR.        &
       rivers_y_bounds_edge(2) > land_y_bounds_edge(2) + dy_toler ) THEN
    CALL log_info( RoutineName,                                                &
                   "There are river gridboxes outside the area of land. " //   &
                   " Runoff to these river sections will be zero." )
  END IF

  ! Look for land gridboxes that have some area outside the river area,
  ! allowing a small tolerance based on land gridbox size.
  dx_toler = frac_toler * land_dx
  dy_toler = frac_toler * land_dy
  IF ( land_x_bounds_edge(1) < rivers_x_bounds_edge(1) - dx_toler  .OR.        &
       land_x_bounds_edge(2) > rivers_x_bounds_edge(2) + dx_toler  .OR.        &
       land_y_bounds_edge(1) < rivers_y_bounds_edge(1) - dy_toler  .OR.        &
       land_y_bounds_edge(2) > rivers_y_bounds_edge(2) + dy_toler ) THEN
    CALL log_info( RoutineName,                                                &
                   "There are land gridboxes outside the area of rivers. " //  &
                   " Runoff from these will not enter the river system." )
  END IF

ELSE

  !----------------------------------------------------------------------------
  ! .NOT. river_regrid
  ! In this case grids are of the same resolution and we have 1:1 mappings.
  !----------------------------------------------------------------------------
  ! Report on how many land and river points match.
  IF ( l_trivial_mapping ) THEN

    CALL log_info( RoutineName,                                                &
                   "Land and river points are identical (trivial mapping)." )

  ELSE

    ! Report on river points for which no matching land point was found.
    num_no_match = COUNT( map_river_to_land_points(:) == 0 )
    IF ( num_no_match == 0 ) THEN
      CALL log_info( RoutineName,                                              &
                     "All river points have a matching land point." )
    ELSE
      CALL log_info( RoutineName,                                              &
                     "Number of river points without a matching land "      // &
                     "point = " // TRIM(to_string(num_no_match))            // &
                     ". Runoff to these river sections will be zero." )
    END IF

    ! Report on land points for which no matching river point was found.
    ! First count river points with matching land points.
    num_match = COUNT( map_river_to_land_points(:) > 0 )
    IF ( num_match == global_land_pts ) THEN
      CALL log_info( RoutineName,                                              &
                     "All land points have a matching river point." )
    ELSE
      CALL log_info( RoutineName,                                              &
                     "Number of land points without a matching river "      // &
                     "point = "                                             // &
                     TRIM(to_string(global_land_pts - num_match))           // &
                     ". Runoff from these will not enter the river system." )
    END IF

  END IF  !  l_trivial_match

END IF  !  rivers_regrid

RETURN
END SUBROUTINE check_river_coverage

!##############################################################################

END MODULE init_rivers_props_mod
#endif
