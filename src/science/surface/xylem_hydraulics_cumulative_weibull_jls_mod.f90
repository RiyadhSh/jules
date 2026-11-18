! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_hydraulics_CW_jls_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName='XYLEM_HYDRAULICS_CW_JLS_MOD'

PUBLIC :: xylem_conductance_CW_jls, leaf_psi_CW_jls

CONTAINS

! *********************************************************************
! Contains routines used to calculate conductance and leaf water
! potential for a cumulative Weibull conductance model.
!
!                  k(psi) = kmax * e^(-(psi / b) ^ c)
!
! *********************************************************************

! ---------------------------------------------------------------------
! Function to calculate xylem conductance using a cumulative weibull
! distribution.
! NOTE: This function is designed for the stomatal optimisation model.
!       It asumes that the input and output arrays are only large
!       enough to contain values for points with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE xylem_conductance_CW_jls( pft,                                      &
                                     n_water_potentials,                       &
                                     open_pnts,                                &
                                     water_potential,                          &
                                  ! INTENT OUT
                                     xylem_conductance                         &
  )

USE pftparm, ONLY:                                                             &
        kmax_pft, conductance_b, conductance_c, kcrit

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, n_water_potentials                                                           &
                            ! Number of water potentials per open point
, open_pnts
                            ! Number of open stomata

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  water_potential(n_water_potentials, open_pnts)
                            ! Water potentials for each open point)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  xylem_conductance(n_water_potentials, open_pnts)
                            ! Xylem conductance for each open point (m/s)

! Local variables
INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_CW_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

xylem_conductance = kmax_pft(pft)                                              &
     * EXP(-ABS(water_potential / conductance_b(pft))                          &
                ** conductance_c(pft))

! Apply minimum conductance limit
xylem_conductance = MAX(xylem_conductance, kcrit(pft))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE xylem_conductance_CW_jls



! ---------------------------------------------------------------------
! Function to calculate the leaf water potential from the transpiration
! rate and rootzone water potential. Uses a cumulative Weibull
! distribution to model conductance.
!
!                  k(psi) = kmax * e^(-(psi / b) ^ c)
!
! The code bellow uses one of two approxiations to estimate the leaf
! water potential from the transpiration rate. This is done to avoid
! the need to intergrate the conductance equation each iteration.
!   1: Applies a zeroth order taylor expansion to the conductance
!      equation (psi_aprox_TE).
!   2: Uses a Newton Raphson aproximation to find the leaf water
!      potential (psi_aprox_NR).
!
! NOTE: This function is designed for the stomatal optimisation model.
!       It returns an output array that only contain values for points
!       with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE leaf_psi_CW_jls( pft,                                               &
                            n_e_leaf,                                          &
                            land_pts,                                          &
                            open_pnts,                                         &
                            veg_index,                                         &
                            open_index,                                        &
                            e_leaf,                                            &
                            root_zone_psi,                                     &
                         ! INTENT OUT
                            leaf_psi,                                          &
                            leaf_k                                             &
  )

USE jules_vegetation_mod, ONLY: som_psi_aprox_method, psi_aprox_TE,            &
                                psi_aprox_NR

USE pftparm, ONLY: kmax_pft, conductance_b, conductance_c, kcrit

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, n_e_leaf                                                                     &
                            ! Number of transpiration values per land point
, land_pts                                                                     &
                            ! Number of land points
, open_pnts                                                                    &
                            ! Number of open stomata
, veg_index(land_pts)                                                          &
                            ! Index of vegetation points on the land grid
, open_index(land_pts)
                            ! Index of open stomata for each land point

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts)                                                  &
                            ! Transpiration rates for each open point (kg/m2/s)
, root_zone_psi(land_pts)
                            ! Water potential in the root zone (Pa)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (m/s)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point (m/s)

! Local variables
INTEGER :: i, j, k, l

! Taylor expansion variables
REAL(KIND=real_jlslsm) :: leaf_k_new
REAL(KIND=real_jlslsm) :: k_conversion_limit
REAL(KIND=real_jlslsm) :: refrence_psi

! Newton-Raphson variables
REAL(KIND=real_jlslsm) :: e_leaf_current(n_e_leaf, open_pnts)
                            ! Current transpiration rates for each open point
                            ! calculated from current leaf_psi. Used in
                            ! Newton-Raphson method. (kg/m2/s)
REAL(KIND=real_jlslsm) :: gamma_root_psi(land_pts)
                            ! Lower incomplete gamma function for the root
                            ! zone water potential. Used in Newton-Raphson
                            ! method.

INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Select the method to approximate the leaf water potential
SELECT CASE (som_psi_aprox_method)

! ---------------------------------------------------------------------
! Apply a zeroth order taylor expansion to simplify the integration
!
! Note: The taylor series expansion is aproximated at the refrence
!       water potential (refrence_psi). The refrence water potential
!       is the mean (midpoint) of the root zone and leaf water
!       potentials.
! ---------------------------------------------------------------------
CASE(psi_aprox_TE)

  ! Calcualte the conductance conversion limit
  k_conversion_limit = 0.1 * kcrit(pft)

  DO j = 1,open_pnts ! Itterate over the open points.
      l = veg_index(open_index(j))

      ! Set the initial leaf water potential to that of the root zone
      leaf_psi(:,j) = root_zone_psi(l)

      ! Calculate the initial xylem conductance for the initial
      !  refrence water potential (this is equal to the root
      !  zone water potential).
      !  k(psi) = kmax * exp(-(psi/b) ** c)
      leaf_k(:,j) = kmax_pft(pft)                                         &
          * EXP( -(leaf_psi(:,j)/conductance_b(pft))**conductance_c(pft) )

      !  Itterate over the transpiration for each openpoint.
      DO i = 1,n_e_leaf

        DO k = 1,100 ! Itterate untill the change in xylem conductance is
                     ! less than k_conversion_limit. or 100 itterations.

          ! Calculate the new leaf water potential.
          leaf_psi(i,j) = root_zone_psi(l) - e_leaf(i,j) / leaf_k(i,j)

          ! Calculate the refrence psi at wich to apply the zeroth order
          ! Taylor series expansion.
          !  psi_ref = mean(leaf_psi, root_zone_psi)
          refrence_psi = 0.5 * (leaf_psi(i,j) + root_zone_psi(l))

          ! Calculate the new xylem conductance
          !  k(psi) = kmax * exp(-(psi/b) ** c)
          leaf_k_new = kmax_pft(pft)                                        &
             * EXP( -(refrence_psi/conductance_b(pft))**conductance_c(pft) )

          ! If the change in xylem conductance is less than
          !  k_conversion_limit then break the loop.
          IF (ABS(leaf_k_new - leaf_k(i,j)) < k_conversion_limit) THEN
            !leaf_k(i,j) = leaf_k_new
            EXIT
          END IF

          ! Set the xylem conductacne to the new value for the next
          ! itteration.
          leaf_k(i,j) = leaf_k_new

        END DO ! Itterate untill conductance converges.
      END DO ! Itterate over the sample ci values.
    END DO ! Itterate over the open points.

    ! Calculate conductance at th leaf water potential. The value calculated
    !  above is the conductance at the refrence psi.
    CALL xylem_conductance_CW_jls(pft, n_e_leaf, open_pnts, leaf_psi, leaf_k)

    leaf_k(:,:) = MAX(leaf_k(:,:), kcrit(pft))

! ---------------------------------------------------------------------
! Apply a Newton-Raphson method to approximate the leaf water potential
! from the transpiration rate. The transpiration rate is given by the
! intigral of the conductance model, k(psi), over the range
! [psi_r, psi_l]. Here the conductance model is given by the cumulative
! Weibull distribution. The integral of k(psi) with respect to psi is;
!
!  E(psi_l, psi_r) = kmax (b/c)
!                     * [ gamma(1/c, (psi_r/b)^c)
!                         - gamma(1/c, (psi_l/b)^c) ]
!
! Where gamma is the lower incomplete gamma function.
!
! The leaf water potential is then calculated by solving the equation;
!
!  0 = E(psi_l, psi_r) - e_leaf
!
! This means that the updated prediction of the leaf water potential
! (psi_l') is given by;
!
!  psi_l' = psi_r - (e_leaf - E(psi_l, psi_r))/k(psi_l)
!
! ---------------------------------------------------------------------
CASE(psi_aprox_NR)

  ! Calculate the lower incomplete gamma function for the root zone water
  ! potential, for all land points.
  gamma_root_psi = incomplete_gamma(land_pts, 1/conductance_c(pft),            &
                                    (root_zone_psi/conductance_b(pft))         &
                                    **conductance_c(pft))

  DO j = 1, open_pnts ! Itterate over the open points.
    l = veg_index(open_index(j))

    ! Do the first iteration outside the loop.

    ! Calculate the xylem conductance at the root zone water potential.
    leaf_k(:,j) = kmax_pft(pft)                                                &
        * EXP( -(root_zone_psi(l)/conductance_b(pft))**conductance_c(pft) )

    ! The intial guess for the leaf water potential is that of the root zone
    ! meaning that the transpiration rate is zero. This simplifies the first
    ! prediction of the leaf water potential to;
    leaf_psi(:,j) = root_zone_psi(l) - e_leaf(:,j) / leaf_k(:,j)

    ! Itterate the water potential
    DO i = 1,2

      ! Calculate the xylem conductance for the current leaf water potential.
      leaf_k(:,j) = kmax_pft(pft)                                              &
          * EXP( -(leaf_psi(:,j)/conductance_b(pft))**conductance_c(pft) )

      ! Calculate the transpiration rate for the current leaf water potential.
      ! First get the incomplete gamma function for the leaf water potential.
      e_leaf_current(:,j) = incomplete_gamma(n_e_leaf,                         &
                                             1/conductance_c(pft),             &
                                             (leaf_psi(:,j)/conductance_b(pft))&
                                              **conductance_c(pft))

      ! Calculate the difference between the incomplete gamma function at the
      ! leaf and root zone water potentials.
      e_leaf_current(:,j) = e_leaf_current(:,j) - gamma_root_psi(l)

      ! Multiply by kmax * (-b/c) to get the transpiration rate
      ! Note: The negative sign is present because conductance_b is negative.
      e_leaf_current(:,j) = e_leaf_current(:,j) * kmax_pft(pft)                &
          * (-conductance_b(pft)/conductance_c(pft))

      ! Calculate the new leaf water potential.
      leaf_psi(:,j) = leaf_psi(:,j)                                            &
             - (e_leaf(:,j) - e_leaf_current(:,j)) / leaf_k(:,j)

    END DO

    ! Calculate the xylem conductance for the current leaf water potential.
    leaf_k(:,j) = kmax_pft(pft)                                              &
          * EXP( -(leaf_psi(:,j)/conductance_b(pft))**conductance_c(pft) )

  END DO

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
        'som_psi_aprox_method should be psi_aprox_TE or psi_aprox_NR')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_psi_CW_jls

! ---------------------------------------------------------------------
! Function to calculate the lower incomplete gamma function. The
! lower incomplete gamma function is defined as;
!
! gamma(a, x) = integral( t^(a-1) * e^(-t) dt, 0, x)
!
! The code here is based on the implimentation in Numerical Recipes
! in FORTRAN 90, page 211, and makes use of the series representation
! of the lower incomplete gamma function.
! ---------------------------------------------------------------------
FUNCTION incomplete_gamma( n_x, a, x ) RESULT( gamma )

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_x ! Size of x array

REAL(KIND=real_jlslsm), INTENT(IN) :: a
REAL(KIND=real_jlslsm), INTENT(IN) :: x(n_x)

REAL(KIND=real_jlslsm) :: gamma(n_x)

! Local variables
REAL(KIND=real_jlslsm) :: step(n_x)
INTEGER :: i

! zhook variables
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='incomplete_gamma'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set initial values
step(:) = 1/a
gamma(:) = step

! Calculate series
DO i = 1,2
  step = step * x / (a + i)
  gamma = gamma + step
END DO

! Multiply by terms outsied the sum
gamma = gamma * EXP(-x) * x**a

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION incomplete_gamma

END MODULE xylem_hydraulics_CW_jls_mod