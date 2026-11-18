! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_hydraulics_SOX_jls_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='XYLEM_HYDRAULICS_SOX_JLS_MOD'

PUBLIC :: xylem_conductance_SOX_jls, leaf_psi_SOX_jls

CONTAINS

! *********************************************************************
! Contains routines used to calculate conductance and leaf water
! potential using the conductance model from SOX, Eller CB et al. 2018.
!
!                                kmax
!                  k(psi) = ---------------
!                            1 + (psi/b)^c
! *********************************************************************

! ---------------------------------------------------------------------
! Function to calculate xylem conductance using the SOX conductance
! model, Eller CB et al. 2018.
! NOTE: This function is designed for the stomatal optimisation model.
!       It asumes that the input and output arrays are only large
!       enough to contain values for points with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE xylem_conductance_SOX_jls( pft,                                    &
                                      n_water_potentials,                     &
                                      open_pnts,                              &
                                      water_potential,                        &
                                   ! INTENT OUT
                                      xylem_conductance                       &
  )

USE pftparm, ONLY:                                                             &
        kmax_pft, conductance_b, conductance_c, kcrit

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                        &
  pft                                                                         &
                            ! Plant functional type index
, n_water_potentials                                                          &
                            ! Number of water potentials per open point
, open_pnts
                            ! Number of open stomata

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  water_potential(n_water_potentials, open_pnts)
                            ! Water potentials for each open point)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  xylem_conductance(n_water_potentials, open_pnts)
                            ! Xylem conductance for each open point (m/s)

! Local variables
INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_SOX_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


xylem_conductance = kmax_pft(pft)                                              &
      / (1 + (water_potential/conductance_b(pft))**conductance_c(pft))


! Apply minimum conductance limit
xylem_conductance = MAX(xylem_conductance, kcrit(pft))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE xylem_conductance_SOX_jls



! ---------------------------------------------------------------------
! Function to calculate the leaf water potential from the transpiration
! rate and rootzone water potential. Uses a cumulative Weibull
! distribution to model conductance.
!
!                                kmax
!                  k(psi) = ---------------
!                            1 + (psi/b)^c
!
! The code bellow uses one of two approxiations to estimate the leaf
! water potential from the transpiration rate. This is done to avoid
! the need to intergrate the conductance equation each iteration.
!   1: Applies a zeroth order taylor expansion to the conductance
!      equation (psi_aprox_TE).
!   2: Uses a precalculated look up table (psi_aprox_NR).
!
! NOTE: This function is designed for the stomatal optimisation model.
!       It returns an output array that only contain values for points
!       with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE leaf_psi_SOX_jls( pft,                                              &
                             n_e_leaf,                                         &
                             land_pts,                                         &
                             open_pnts,                                        &
                             veg_index,                                        &
                             open_index,                                       &
                             e_leaf,                                           &
                             root_zone_psi,                                    &
                          ! INTENT OUT
                             leaf_psi,                                         &
                             leaf_k                                            &
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

! Taylor series variables
REAL(KIND=real_jlslsm) :: leaf_k_new
REAL(KIND=real_jlslsm) :: k_conversion_limit
REAL(KIND=real_jlslsm) :: refrence_psi

! Newton-Raphson variables
REAL(KIND=real_jlslsm) :: e_leaf_current(n_e_leaf, open_pnts)
                            ! Current transpiration rates for each open point
                            ! calculated from current leaf_psi. Used in
                            ! Newton-Raphson method. (kg/m2/s)
REAL(KIND=real_jlslsm) :: root_psi_component(land_pts)
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
      leaf_k(:,j) = kmax_pft(pft)                                              &
                 / (1 + (leaf_psi(:,j)/conductance_b(pft))**conductance_c(pft))

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
          leaf_k_new = kmax_pft(pft)                                           &
                 / (1 + (refrence_psi/conductance_b(pft))**conductance_c(pft))

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
    CALL xylem_conductance_SOX_jls(pft, n_e_leaf, open_pnts, leaf_psi, leaf_k)

    leaf_k(:,:) = MAX(leaf_k(:,:), kcrit(pft))

! ---------------------------------------------------------------------
! Apply a Newton-Raphson method to approximate the leaf water potential
! from the transpiration rate. The transpiration rate is given by the
! intigral of the conductance model, k(psi), over the range
! [psi_r, psi_l]. Here the conductance model is that used in the SOX
! 2018 paper. The integral of k(psi) with respect to psi is;
!
!  E(psi_l, psi_r) = kmax (b/c)
!                     * [ psi_l * 2F1(1, 1/c; 1 + 1/c; -(psi_l/b)^c)
!                        - psi_r * 2F1(1, 1/c; 1 + 1/c; -(psi_r/b)^c)]
!
! Where 2F1() is the Gaussian Hypergeometric series. (This integral
! was calculated using Wolfram Alpha's online integrator)
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

  ! Calculate the lower 2F1 function for the root zone water
  ! potential, for all land points.
  root_psi_component = SOX_2F1(land_pts,                                       &
                               1 + 1/conductance_c(pft),                       &
                -(root_zone_psi/conductance_b(pft))**conductance_c(pft))

  root_psi_component = root_zone_psi * root_psi_component

  DO j = 1, open_pnts ! Itterate over the open points.
    l = veg_index(open_index(j))

    ! Do the first iteration outside the loop. By setting the initial guess
    ! to the root zone water potential, the transpiration rate is zero.

    ! Calculate the xylem conductance at the root zone water potential.
    leaf_k(:,j) = kmax_pft(pft)                                                &
        / (1 + (root_zone_psi(l)/conductance_b(pft))**conductance_c(pft))

    ! The intial guess for the leaf water potential is that of the root zone
    ! meaning that the transpiration rate is zero. This simplifies the first
    ! prediction of the leaf water potential to;
    leaf_psi(:,j) = root_zone_psi(l) - e_leaf(:,j) / leaf_k(:,j)

    !WRITE(*,*) ' '
    !WRITE(*,*) 'root_psi', root_zone_psi(l)
    !WRITE(*,*) 'leaf_psi', leaf_psi(400,j)

    ! Itterate the water potential
    DO i = 1,2

      ! Calculate the xylem conductance for the current leaf water potential.
      leaf_k(:,j) = kmax_pft(pft)                                              &
        / (1 + (leaf_psi(:,j)/conductance_b(pft))**conductance_c(pft))

      ! Calculate the transpiration rate for the current leaf water potential.
      ! First get the incomplete gamma function for the leaf water potential.
      e_leaf_current(:,j) = SOX_2F1(n_e_leaf,                                  &
                                    1 + 1/conductance_c(pft),                  &
                   -(leaf_psi(:,j)/conductance_b(pft))**conductance_c(pft))

      e_leaf_current(:,j) = leaf_psi(:,j) * e_leaf_current(:,j)

      ! Calculate the difference between the incomplete gamma function at the
      ! leaf and root zone water potentials.
      e_leaf_current(:,j) = e_leaf_current(:,j) - root_psi_component(l)

      ! Multiply by kmax to get the transpiration rate.
      e_leaf_current(:,j) = - e_leaf_current(:,j) * kmax_pft(pft)

      ! Calculate the new leaf water potential.
      leaf_psi(:,j) = leaf_psi(:,j)                                            &
             + (e_leaf(:,j) - e_leaf_current(:,j)) / leaf_k(:,j)

      !WRITE(*,*) 'leaf_psi', leaf_psi(400,j)

    END DO

    ! Calculate the xylem conductance for the current leaf water potential.
    leaf_k(:,j) = kmax_pft(pft)                                                &
        / (1 + (leaf_psi(:,j)/conductance_b(pft))**conductance_c(pft))

  END DO

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
        'som_psi_aprox_method should be psi_aprox_TE or psi_aprox_NR')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_psi_SOX_jls

! ---------------------------------------------------------------------
! Function to calculate the Gaussian Hypergeometric series
! 2F1(af,bf;cf;xf). THis function is used to calculate the
! transpiration rate from the leaf water potential for the SOX
! conductance model.
!
! For the SOX conductance model the inputs to the 2F1 function are
! allways;
!   af = 1
!   bf = 1/c
!   cf = 1 + 1/c
!   xf = -(psi/b)^c <= 0
!
! The series solution for 2F1(af,bf;cf;xf) is;
!
!                            (af)_n * (bf)_n     xf^n
!   2F1(af,bf;cf;xf) = SUM[ ----------------- * ------ , n=0,inf]
!                                (cf)_n           n!
!
! over the range |xf| < 1. The notation (q)_n is the Pochhammer symbol,
! usd to denote the rising factorial;
!
!    (q)_n = q * (q+1) * (q+2) * ... * (q+n-1)
!
!    (q)_0 = 1
!
! To allow for the clculation of 2F1 for all xf <=0 the code below uses
! the analytic continuation of 2F1 described in equation B.5 of
! W. Becken, P. Schmelcher 1997 (The analytic continuation of the
! Gaussian hypergeometric function 2F1(a,b;c;z) for arbitrary
! parameters). Equation B.5 takes the form;
!
! 2F1(af,bf;cf;xf) / Gamma(cf) =
!      (1 - xf)^(-af) * 2F1(af,cf-bf;cf;xf/(xf-1))/Gamma(cf)
!
! Note that this equation can be used because cf - af - bf = 0.
! Canceling out the two Gamma functions and substituting for the values
! of, bf and cf.
!
! 2F1(1,1/c;1+1/c;xf) =
!     (1 + xf)^(-1) * 2F1(1,1;1+1/c;xf/(xf - 1))
!
! Given that xf <= 0, xf/(xf-1) exists over the range [0,1). Meaning
! that the series solution for 2F1(1,1;1+1/c;xf/(xf-1)) can be
! calculatedusing the series solution for 2F1.
!
!  2F1(1,1;1+1/c;xf/(xf-1)) =
!
!          (1)_n * (1)_n       xf         1
!    SUM[ --------------- * (------)^n * ---- , n=0,inf]
!            (1+1/c)_n        xf-1        n!
!
! Since (1)_n = n! this simplifies to;
!
!  2F1(1,1;1+1/c;xf/(xf-1)) =
!
!              n!          xf
!    SUM[ ----------- * (------)^n, n=0,inf]
!          (1+1/c)_n      xf-1
!
! The value within this sum would be expensive to calculate for large
! n. However, it is possible to re-write the value within the sum such
! that it is a multiple ofthe previous itteration (n-1);
!
!  2F1(1,1;1+1/c;xf/(xf-1)) =
!
!                       xf          (n-1)!           xf
!   1 + SUM[ n * c * (------) * --------------- * (------)^(n-1), n=1,inf]
!                      xf-1      (1+1/c)_(n-1)      xf-1
!
! Note the first term in the sum (n=0) is equal to one and has been
! moved outside the sum.
!
! If we let Sn be the nth term in the sum, with S0 = 1, then the nth
! term in the sum is given by;
!
!                   xf
!   Sn = n * c * (------) * Sn-1
!                  xf-1
!
! Where Sn-1 is the previous term in the sum.
! ---------------------------------------------------------------------
FUNCTION SOX_2F1(n_xf, cf, xf) RESULT( sum )

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_xf ! Size of xf array

REAL(KIND=real_jlslsm), INTENT(IN) :: cf
REAL(KIND=real_jlslsm), INTENT(IN) :: xf(n_xf)

REAL(KIND=real_jlslsm) :: sum(n_xf)

! Local variables
REAL(KIND=real_jlslsm) :: current_term(n_xf)
INTEGER :: i

! zhook variables
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='SOX_2F1'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the 0th term in the sum.
current_term = 1.0
sum = 1.0

! Itterate over the sum
DO i = 1, 2
  current_term = current_term * i * cf * xf / (xf - 1)
  sum = sum + current_term
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION SOX_2F1

END MODULE xylem_hydraulics_SOX_jls_mod