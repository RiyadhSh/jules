! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE stom_opt_jls_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='STOM_OPT_JLS_MOD'

PRIVATE stom_opt_mod_ci
PUBLIC stom_opt_mod

CONTAINS

! *********************************************************************
! Contains routines to calculate the optimal stomatal conductance and
! net carbon uptake.
!
! Refrences:
! TODO: Add refrence to profit max paper.
! TODO: ADD refrence to JULES paper.
! *********************************************************************

SUBROUTINE stom_opt_mod (                                                      &
! IN
        land_pts, som_base_parm, pft, open_pts, open_index,                    &
        pft_photo_model, veg_index,                                            &
        ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,           &
        km, dq, qs, je, t_leaf,                                                &
! IN OUT
        rd,                                                                    &
! OUT
        ci, al, el, flux_o3, fo3, gl, psi_leaf,                                &
        carbon_gain_out, hydraulic_cost_out, leaf_k                            &
)

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE jules_vegetation_mod, ONLY:                                                &
        som_base_parm_ci, som_base_parm_psi, som_n_sample,                     &
        profit_max_profit_model, SOX_profit_model, som_profit_model

USE pftparm, ONLY:                                                             &
        min_gl_pft, kcrit, kcrit_fractional_loss, kmax_pft

USE xylem_hydraulics_jls_mod, ONLY: xylem_conductance_jls

!-----------------------------------------------------------------------------
! IN integer variables.
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  land_pts                                                                     &
!                           ! Number of land points to be processed
, som_base_parm                                                                &
!                           ! Switch used to determin the physical property
!                           ! that is uniformly itterated over when
!                           ! determaning the optimal stomatal conductance.
!                           !     1 => leaf intercellular carbon (Ci)
!                           !     2 => leaf water potential (psi)
, pft                                                                          &
!                           ! Plant functional type
,open_pts                                                                      &
                            ! Number of land points with open stomata.
,open_index(land_pts)                                                          &
                            ! Index of land points with open stomata.
,pft_photo_model                                                               &
                            ! Indicates which photosynthesis model to use for
                            ! the current PFT.
,veg_index(land_pts)
                            ! Index of vegetated points on the land grid.

!-----------------------------------------------------------------------------
! IN real variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 ca(land_pts)                                                                  &
                             ! Atmostpheric CO2 pressure (Pa)
,psi_root_zone(land_pts)                                                       &
                             ! Root zone water potential (Pa)
,acr(land_pts)                                                                 &
                             ! Absorbed PAR (mol photons/m2/s).
,apar(land_pts)                                                                &
                             ! Absorbed PAR (W m-2).
,oi(land_pts)                                                                  &
                             ! Intercellular oxygen pressure (Pa).
,vcmax(land_pts)                                                               &
                             ! Maximum rate of carboxylation of Rubisco
                             !  (mol CO2/m2/s).
,kc(land_pts)                                                                  &
                             ! Michaelis-Menten constant for CO2 (Pa).
,ko(land_pts)                                                                  &
                             ! Michaelis-Menten constant for O (Pa).
,ccp(land_pts)                                                                 &
                            ! Photorespiratory compensatory point (Pa).
,pstar(land_pts)                                                               &
                            ! Atmospheric pressure (Pa).
,km(land_pts)                                                                  &
                            ! A combination of Michaelis-Menten and other
                            ! terms. Only used with the Farquhar model.
,dq(land_pts)                                                                  &
                            ! Specific humidity deficit of air (kg H2O/kg air).
                            ! Leaf saturation vapour pressure (Pa).
,qs(land_pts)                                                                  &
                            ! Saturated specific humidity
                            !      (kg H2O/kg air).
,je(land_pts)                                                                  &
                            ! Electron transport rate (mol m-2 s-1)
,t_leaf(land_pts)
                            ! Leaf temperature (K)

!-----------------------------------------------------------------------------
! IN OUT real variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 rd(land_pts)              ! Dark respiration (mol CO2/m2/s).
                             ! This is modified only if l_o3_damage=.TRUE..

!-----------------------------------------------------------------------------
! OUT real variables.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 ci(land_pts)                                                                  &
                            ! Internal CO2 pressure (Pa).
,al(land_pts)                                                                  &
                            ! Net Leaf photosynthesis (mol CO2/m2/s).
,el(land_pts)                                                                  &
                            ! Transpiration rate (mol H2O/m2/s)
,flux_o3(land_pts)                                                             &
                            ! Flux of O3 to stomata (nmol O3/m2/s).
,fo3(land_pts)                                                                 &
                            ! Ozone exposure factor.
,gl(land_pts)                                                                  &
                            ! Leaf conductance for H2O (m/s).
,psi_leaf(land_pts)                                                            &
                            ! Leaf water potential (Pa)
,leaf_k(land_pts)
                            ! Xylem conductance at leaf water potential (m/s)

! TEMPORARY: output variables for testing
REAL(KIND=real_jlslsm) ::                                                      &
 carbon_gain_out(land_pts)                                                     &
                            ! Carbon gain for each leaf state
,hydraulic_cost_out(land_pts)
                            ! Hydraulic cost for each leaf state

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
 optimal_index                                                                 &
                            ! Holds index of the optimal stomatal conductance
                            !  for each land point.
,i,j,l                                                                         &
                            ! Itterators
,errcode
                            ! Error code to pass to ereport.

LOGICAL ::                                                                     &
 l_good_sample(som_n_sample, open_pts)
                            ! Logical used to determin if the current sample
                            !  is physicaly reasonable and so can be "chosen"
                            !  by the leaf.

! Arrays containing results of sampling over the base parameter. Only need to
!  do so for locations with open stomata
REAL(KIND=real_jlslsm) ::                                                      &
 ci_sample(0:som_n_sample, open_pts)                                           &
                            ! Internal CO2 pressure (Pa).
,al_sample(0:som_n_sample, open_pts)                                           &
                            ! Net Leaf photosynthesis (mol CO2/m2/s).
,gl_sample(0:som_n_sample, open_pts)                                           &
                            ! Leaf conductance for H2O (m/s).
,kl_sample(0:som_n_sample, open_pts)                                           &
                            ! Xylem conductance at leaf (m/s).
,kl_SOX(0:som_n_sample, open_pts)                                              &
                            ! Xylem conductance for SOX profit model.
                            !  Calculated from the average of the leaf and
                            !  root zone water potentials. (m/s)
,psi_sample(0:som_n_sample, open_pts)                                          &
                            ! Water potential at leaf (Pa)
,SOX_mean_psi(0:som_n_sample, open_pts)                                        &
                            ! Mean water potential between leaf and root zone
                            ! (Pa). Used in SOX profit model.
,el_sample(0:som_n_sample, open_pts)                                           &
                            ! Transpiration rate (mol H2O/m2/s)
,profit(0:som_n_sample, open_pts)                                              &
                            ! Profit for each leaf state
,carbon_gain(0:som_n_sample, open_pts)                                         &
                            ! Carbon gain for each leaf state
,hydraulic_cost(0:som_n_sample, open_pts)                                      &
                            ! Hydraulic cost for each leaf state
,max_al(open_pts)                                                              &
                            ! Maximum net leaf photosynthesis (mol CO2/m2/s)
,max_kl(open_pts)
                            ! Maximum xylem conductance (m/s)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='STOM_OPT_MOD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------------
!  Parameter setup
! ----------------------------------------------------------------------------
ci(:)      = ca(:)
! Set to -rd to get net photosynthesis for points with closed stomata. This
! is overwriten for points with open stomata.
al(:)      = -rd(:)
el(:)      = 0.0
flux_o3(:) = 0.0
fo3(:)     = 0.0
gl(:)      = min_gl_pft(pft)
psi_leaf(:)= psi_root_zone(:)
l_good_sample(:,:) = .TRUE.
leaf_k(:)  = kmax_pft(pft)

carbon_gain(:,:) = 0.0
hydraulic_cost(:,:) = 0.0
profit(:,:) = 0.0

! If there are no land points with open stomata then no calculation is needed.
IF(0 == open_pts) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  return
END IF

! ----------------------------------------------------------------------------
!  Set input environmental properties for testing
!  TODO: Remove this section
! ----------------------------------------------------------------------------

! Set temporary output variables to unrealistic values.
! Note hydraulic cost is set to 2 such that both the value is obveously wrong
! and the profit calculated would be negative.
carbon_gain_out(:) = -1.0
hydraulic_cost_out(:) = 2.0

! ----------------------------------------------------------------------------
!  Calculate leaf properties as a function of the selected paramiter,
!    eather ci or psi.
! ----------------------------------------------------------------------------
SELECT CASE ( som_base_parm )
  CASE (som_base_parm_ci)
    CALL stom_opt_mod_ci(                                                      &
    ! IN
        land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,       &
        rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,       &
        km, dq, qs, je, t_leaf,                                                &
    ! OUT
        ci_sample, al_sample, gl_sample, kl_sample, psi_sample, el_sample      &
            )

  CASE (som_base_parm_psi)
    WRITE(*,*) "psi not implimented yet"

  CASE DEFAULT
    errcode = 101  !  a hard error
    CALL ereport(RoutineName, errcode,                                         &
               'som_base_parm shoud be 1 or 2')
END SELECT

! ----------------------------------------------------------------------------
!  Make mask of those samples that are physicaly reasonable.
! ----------------------------------------------------------------------------
!TODO: add open mp to loop
DO j = 1, open_pts
  l = open_index(j)

  ! Requier psi <= psi_root_zone
  l_good_sample(:,j) = l_good_sample(:,j)                                    &
          .and. psi_sample(1:,j) <= psi_root_zone(l) + TINY(psi_root_zone(l))

  ! Requier k > k_crit
  l_good_sample(:,j) = l_good_sample(:,j)                                    &
                       .and. kl_sample(1:,j) > kcrit(pft)

  ! Requier ci < ca
  !l_good_sample(i,j) = l_good_sample(i,j) .and. ci_sample(i,j) < ca(l)
END DO

! ----------------------------------------------------------------------------
!  Calculate the profit for each leaf state.
!    Profit models:
!      1: profit_max_profit_model
!           CG = An/max(An)
!           HC = (ki_max - k)/(ki_max - k_crit)
!           Profit = CG-HC
!
!      2: SOX_profit_model
!           CG = An
!           HC = 1 - (k - k_crit)/(ki_max - k_crit)
!           Profit = CG*(1-HC)
! ----------------------------------------------------------------------------
!TODO: add open mp to loop
SELECT CASE (som_profit_model)
  CASE (profit_max_profit_model)

    ! Get the maximum net leaf photosynthesis and xylem conductance for each land
    ! point.
    max_al = MAXVAL(al_sample(1:,:), MASK = l_good_sample, DIM = 1)
    max_kl = MAXVAL(kl_sample(1:,:), MASK = l_good_sample, DIM = 1)

    ! Loop over land points with open stomata
    DO j = 1, open_pts

      ! CG = An/max(An)
      carbon_gain(:,j) = al_sample(:,j)/max_al(j)

      ! HC = (ki_max - k)/(ki_max - k_crit)
      hydraulic_cost(:,j) =   (max_kl(j) - kl_sample(:,j))                     &
                            / (max_kl(j) - kcrit(pft))

      ! Profit = CG-HC
      profit(:,j) = carbon_gain(:,j) - hydraulic_cost(:,j)

    END DO

  CASE (SOX_profit_model)

    ! Sox uses the xylem conductance for the average of the
    !  leaf and root zone water potentials.
    DO j = 1, open_pts
      l = open_index(j)

      SOX_mean_psi(:,j) = 0.5 * (psi_sample(:,j) + psi_root_zone(l))
    END DO

    CALL xylem_conductance_jls(pft, som_n_sample, open_pts,                    &
                               SOX_mean_psi,                                   &
                               kl_SOX)

    ! Get the maximum xylem conductance for each land point.
    max_kl = MAXVAL(kl_SOX(1:,:), MASK = l_good_sample, DIM = 1)

    ! CG = An
    carbon_gain(:,:) = al_sample(:,:)

    ! Loop over land points with open stomata
    DO j = 1, open_pts
        ! HC = 1 - k/ki_max
        hydraulic_cost(:,j) = 1 - (kl_SOX(:,j)-kcrit(pft))                     &
                                  /(max_kl(j)-kcrit(pft))

        ! Profit = CG*(1-HC)
        profit(:,j) = carbon_gain(:,j)*(1-hydraulic_cost(:,j))
    END DO ! i Land points

  CASE DEFAULT
    errcode = 101  !  a hard error
    CALL ereport(RoutineName, errcode,                                         &
               'profit_model shoud be profit_max or SOX')

END SELECT ! profit_model

! ----------------------------------------------------------------------------
!  Identify maximum profit and set leaf properties to those of the optimal
!   state.
!
!  NOTE: The zeroth entry of the sample arrays corresponds to closed stomata.
!        i.e. ci = ca, al = 0, gl = 0, kl = 0, psi = psi_root_zone, el = 0.
!        If the mask (l_good_sample) is false for all samples then MAXLOC
!        will return 0 for the optimal index. The leaf properties will then
!        be set to those of the closed stomata state.
! ----------------------------------------------------------------------------
DO j = 1, open_pts
  l = open_index(j)

  !Find the index of the maximum profit, excluding the closed stomata state.
  ! Returns 0 if all samples are not physicaly reasonable.
  optimal_index = MAXLOC(profit(1:,j), DIM = 1, MASK = l_good_sample(:,j))

  !Set the leaf properties to those of the optimal state
  ci(l) = MAX(0.0, ci_sample(optimal_index,j))
  al(l) = MAX(0.0, al_sample(optimal_index,j))
  gl(l) = MAX(0.0, gl_sample(optimal_index,j))
  psi_leaf(l) = psi_sample(optimal_index,j)
  el(l) = el_sample(optimal_index,j)
  leaf_k(l) = kl_sample(optimal_index,j)

  ! Set the output variables for testing
  carbon_gain_out(l) = carbon_gain(optimal_index,j)
  hydraulic_cost_out(l) = hydraulic_cost(optimal_index,j)
END DO


!TODO: Redetermin if stomata are open or closed?

!TODO: Calculate ozone exposure here. Create subroutine using code in
!TODO:  leaf_jls_mod.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE stom_opt_mod

!-----------------------------------------------------------------------------
! Stomatal optimisation model that iterates over leaf ci
!-----------------------------------------------------------------------------
SUBROUTINE stom_opt_mod_ci(                                                    &
! IN
        land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,       &
        rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp,              &
        pstar, km, dq, qs, je, t_leaf,                                         &
! OUT
        ci_sample, al_sample, gl_sample, kl_sample, psi_sample,el_sample       &
)

USE xylem_hydraulics_jls_mod, ONLY: leaf_psi_jls

USE jules_vegetation_mod, ONLY:                                                &
        som_n_sample, photo_collatz, photo_farquhar, CW_conductance,           &
        SOX_conductance, som_psi_aprox_method

USE pftparm, ONLY:                                                             &
        leaf_crit, c3, alpha, pft_conductance_model,                           &
        kmax_pft, conductance_b, conductance_c, kcrit

USE jules_surface_mod, ONLY: fwe_c3, fwe_c4
USE jules_surface_mod, ONLY: beta1, beta2, ratio, ratio_o3

USE planet_constants_mod, ONLY: repsilon, one_minus_epsilon

USE c_rmol, ONLY: rmol

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

!-----------------------------------------------------------------------------
! IN integer variables
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  land_pts                                                                     &
!                           ! Number of land points to be processed
, pft                                                                          &
                            ! Plant functional type
, open_pts                                                                     &
                            ! number of land points with open stomata
, open_index(land_pts)                                                         &
                            ! Index of each land point with open stmata
, pft_photo_model                                                              &
                            ! Indicates which photosynthesis model to use for
                            ! the current PFT.
, veg_index(land_pts)
                            ! Index of vegetated points on the land grid.

!-----------------------------------------------------------------------------
! IN real variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 rd(land_pts)                                                                  &
                             ! Dark respiration (mol CO2/m2/s).
,ca(land_pts)                                                                  &
                             ! Atmostpheric CO2 pressure (Pa)
,psi_root_zone(land_pts)                                                       &
                             ! Water potential in the root zone (Pa)
,acr(land_pts)                                                                 &
                             ! Absorbed PAR (mol photons/m2/s).
,apar(land_pts)                                                                &
                             ! Absorbed PAR (W m-2).
,oi(land_pts)                                                                  &
                             ! Intercellular oxygen pressure (Pa).
,vcmax(land_pts)                                                               &
                             ! Maximum rate of carboxylation of Rubisco
                             !  (mol CO2/m2/s).
,kc(land_pts)                                                                  &
                             ! Michaelis-Menten constant for CO2 (Pa).
,ko(land_pts)                                                                  &
                             ! Michaelis-Menten constant for O (Pa).
,ccp(land_pts)                                                                 &
                             ! Photorespiratory compensatory point (Pa).
,pstar(land_pts)                                                               &
                             ! Atmospheric pressure (Pa).
,km(land_pts)                                                                  &
                             ! A combination of Michaelis-Menten and other
                             ! terms. Only used with the Farquhar model.
,dq(land_pts)                                                                  &
                            ! Specific humidity deficit of air (kg H2O/kg air).
                            ! Leaf saturation vapour pressure (Pa).
,qs(land_pts)                                                                  &

                            ! Saturated specific humidity
                            !      (kg H2O/kg air).
,je(land_pts)                                                                  &
                             ! Electron transport rate (mol m-2 s-1)
,t_leaf(land_pts)
                             ! Leaf temperature (K)

!-----------------------------------------------------------------------------
! OUT real variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 ci_sample(0:som_n_sample, open_pts)                                             &
                            ! Internal CO2 pressure (Pa).
,al_sample(0:som_n_sample, open_pts)                                             &
                            ! Net Leaf photosynthesis (mol CO2/m2/s).
,gl_sample(0:som_n_sample, open_pts)                                             &
                            ! Leaf conductance for H2O (m/s).
,kl_sample(0:som_n_sample, open_pts)                                             &
                            ! Xylem conductance at leaf (m/s).
,psi_sample(0:som_n_sample, open_pts)                                            &
                            ! Water potential at leaf (Pa)
,el_sample(0:som_n_sample, open_pts)
                            ! Transpiration rate (mol H2O/m2/s)

!-----------------------------------------------------------------------------
! Local integer variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
 i,j,k,l                                                                       &
                            ! Internal itterators
,errcode
                            ! Error code to pass to ereport.

!-----------------------------------------------------------------------------
! Local real variables.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 wcarb_sample(0:som_n_sample, open_pts)                                          &
                            ! Rubisco limited photosynthesis
!                           ! rate (mol CO2/m2/s).
,wlite_sample(0:som_n_sample, open_pts)                                          &
                            ! Electron transport limited photosynthesis
!                           ! rate (mol CO2/m2/s).
,wexpt_sample(0:som_n_sample, open_pts)                                          &
                            ! Export-limited gross photosynthetic rate
!                           ! (mol CO2/m2/s). Not used with Farquhar model.
,wl_sample(0:som_n_sample, open_pts)                                             &
                            ! Gross photosynthetic rate (mol CO2/m2/s).
,b1(0:som_n_sample, open_pts)                                                    &
,b2(0:som_n_sample, open_pts)                                                    &
,b3(0:som_n_sample, open_pts)                                                    &
                            ! Temporary variables for use in the quadratic
                            ! equation.
,beta1p2m4, beta2p2m4                                                          &
                            ! beta[12] ** 2 * 4.
,kl_sample_new, psi_sample_new                                                 &
                            ! Temporary variables for use when estimating
                            ! leaf water potential from transpiration rate.
,vpd(land_pts)                                                                 &
                            ! Vapour pressure deficit (Pa)
,vpd_factor                                                                    &
                            ! Factor used to convert specific humidity to
                            ! vapour pressure.
,glco2_sample(0:som_n_sample, open_pts)                                        &
                            ! Leaf conductance for CO2 (m/s)
,conv(land_pts)                                                                &
                            ! Factor for converting mol/m3 into Pa (J/m3)
,conductance_conversion_limit
                            ! Limit for the change in conductance when
                            ! estimating leaf water potential from
                            ! transpiration rate.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='STOM_OPT_MOD_CI'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Setup parameters.
! NOTE: The zeroth entry to the sample array for each land point corresponds
!        to closed stomata and is only used when determaning the optimal leaf
!        state in stom_opt_mod above. This is why all loops are over the
!        range 1:som_n_sample, NOT 0:som_n_sample.
!-----------------------------------------------------------------------------

al_sample(:,:)      = 0.0
gl_sample(:,:)      = 0.0
kl_sample(:,:)      = 0.0
el_sample(:,:)      = 0.0

DO i = 1, som_n_sample
    !NOTE: Array entry 1 is at the CO2 compensation point (hence i-1), and
    !       the last step is less than the CO2 concentration point in the
    !       air to avoid infinite transpiration rates.
    ci_sample(i,:) = MAX(ccp(:),0.0)                                &
                     + (i-1) * (ca(:) - MAX(ccp(:),0.0))/som_n_sample
  END DO
!-----------------------------------------------------------------------------
! Calculate photosynthetic rates.
! Code based on that in leaf_limits_mod.F90 and leaf.F90 with loops over
! sample_ci added.
!-----------------------------------------------------------------------------
SELECT CASE ( pft_photo_model )
CASE ( photo_collatz )
  !---------------------------------------------------------------------------
  ! Use the Collatz models (for C3 or C4 plants).
  !---------------------------------------------------------------------------

  IF (c3(pft) == 1) THEN

!$OMP PARALLEL DO IF(open_pts > 1)                                             &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,j,i)                                                           &
!$OMP SHARED(open_pts,veg_index,open_index,wcarb_sample,vcmax,ci_sample,       &
!$OMP        ccp,kc,oi,ko,wlite_sample,pft,wexpt_sample,fwe_c3,alpha,acr,      &
!$OMP        som_n_sample)
    DO j = 1,open_pts
      l = veg_index(open_index(j))
      ! The numbers in these equations are from Cox, HCTN 24,
      ! "Description ... Vegetation Model", equations 54 and 55.
      wcarb_sample(:,j) = vcmax(l) * (ci_sample(:,j) - ccp(l))                 &
                          / (ci_sample(:,j) + kc(l) * (1.0 + oi(l) / ko(l)))
      wlite_sample(:,j) = alpha(pft) * acr(l)                                  &
                          * (ci_sample(:,j) - ccp(l))                          &
                          / (ci_sample(:,j) + 2.0 * ccp(l))
      wlite_sample(:,j) = MAX(wlite_sample(:,j), TINY(1.0e0))
      wexpt_sample(:,j) = fwe_c3 * vcmax(l)
    END DO
!$OMP END PARALLEL DO

  ELSE
    !  C4
!$OMP PARALLEL DO IF(open_pts > 1)                                             &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,j,i)                                                           &
!$OMP SHARED(open_pts,veg_index,open_index,wcarb_sample,vcmax,wlite_sample,    &
!$OMP        pft,wexpt_sample,pstar,alpha,fwe_c4,ci_sample,acr,som_n_sample)
    DO j = 1,open_pts
      l = veg_index(open_index(j))
      wcarb_sample(:,j) = vcmax(l)
      wlite_sample(:,j) = alpha(pft) * acr(l)
      wlite_sample(:,j) = MAX(wlite_sample(:,j), TINY(1.0e0))
      wexpt_sample(:,j) = fwe_c4 * vcmax(l) * ci_sample(:,j) / pstar(l)
    END DO
!$OMP END PARALLEL DO

  END IF  !  c3

  !---------------------------------------------------------------------------
  ! Calculate the co-limited rate of gross photosynthesis.
  !---------------------------------------------------------------------------
  beta1p2m4 = 4 * beta1 * beta1
  beta2p2m4 = 4 * beta2 * beta2


  b1(:,:) = beta1
  b2(:,:) = - (wcarb_sample(:,:) + wlite_sample(:,:))
  b3(:,:) = wcarb_sample(:,:) * wlite_sample(:,:)

  wl_sample(:,:) = -b2(:,:) / (2.0 * b1(:,:))                                  &
          - SQRT(b2(:,:) * b2(:,:) / beta1p2m4      - b3(:,:) / b1(:,:))

  b1(:,:) = beta2
  b2(:,:) = - (wl_sample(:,:) + wexpt_sample(:,:))
  b3(:,:) = wl_sample(:,:) * wexpt_sample(:,:)

  wl_sample(:,:) = -b2(:,:) / (2.0 * b1(:,:))                                  &
          - SQRT(b2(:,:) * b2(:,:) / beta2p2m4       - b3(:,:) / b1(:,:))


CASE ( photo_farquhar )
  !---------------------------------------------------------------------------
  ! Use the Farquhar model (for C3 plants only).
  !---------------------------------------------------------------------------
!$OMP PARALLEL DO IF(open_pts > 1)                                             &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,j,i)                                                           &
!$OMP SHARED(open_pts,veg_index,open_index,wcarb_sample,vcmax,wlite_sample,    &
!$OMP        ci_sample,ccp,km,je,som_n_sample)
  DO j = 1,open_pts
    l = veg_index(open_index(j))
      wcarb_sample(:,j) = vcmax(l) * ( ci_sample(:,j) - ccp(l) )               &
                          / ( ci_sample(:,j) + km(l) )
      wlite_sample(:,j) = je(l) / 4.0 * ( ci_sample(:,j) - ccp(l) )            &
                          / ( ci_sample(:,j) + 2.0 * ccp(l) )
      wlite_sample(:,j) = MAX(wlite_sample(:,j), TINY(1.0e0))
  END DO
!$OMP END PARALLEL DO

  !---------------------------------------------------------------------------
  ! Calculate minimum of Rubisco and electron transport limited rates.
  !---------------------------------------------------------------------------
  wl_sample(:,:) = MIN( wcarb_sample(:,:), wlite_sample(:,:) )


CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
               'pft_photo_model should be photo_collatz or photo_farquhar')

END SELECT  !  pft_photo_model

! ----------------------------------------------------------------------------
! Calculate net photosynthetic rate, mol CO2 /m2/s.
!  al = wl - rd
! ----------------------------------------------------------------------------
!TODO: add open mp to loop
DO j = 1,open_pts
  l = veg_index(open_index(j))

  al_sample(:,j) = wl_sample(:,j) - rd(l)
END DO

! ----------------------------------------------------------------------------
! Calculate stomatal conductance in units of m/s.
!  gsw = (1/1.6)*(A/(ca-ci))
! ----------------------------------------------------------------------------
!TODO: add open mp to loop
DO j = 1,open_pts
  l = veg_index(open_index(j))

  !---------------------------------------------------------------------------
  ! Calculate the factor for converting mol/m3 into Pa (J/m3).
  !---------------------------------------------------------------------------
  conv(l) = rmol * t_leaf(l)

  !---------------------------------------------------------------------------
  ! Diagnose the leaf conductance
  !---------------------------------------------------------------------------
  glco2_sample(:,j) = (al_sample(:,j) * conv(l)) / (ca(l) - ci_sample(:,j))
  gl_sample(:,j)    = ratio * glco2_sample(:,j)

END DO

! ----------------------------------------------------------------------------
! Calculate transpiration rate.
!   E = vpd * gl / pstar
! ----------------------------------------------------------------------------
!TODO: add open mp to loop

! Calculate the factor for converting specific humidity to vapour pressure.
! Taken from leaf_limits_mod.F90, the 1.0e3 has been removed to give vpd
! in Pa.
vpd_factor = 1.0 / repsilon

DO j = 1,open_pts
  l = veg_index(open_index(j))

  ! Convert the specific humidity to vapour pressure. This formula is taken
  ! from leaf_limits_mod.F90.
  vpd(l) = dq(l) * pstar(l) * vpd_factor

  ! Calcualte transpiration rate in m/s
  el_sample(:,j) = vpd(l) * gl_sample(:,j) / pstar(l)

  ! Convert from m/s to mol H2O/m2/s
  el_sample(:,j) = el_sample(:,j) * pstar(l) / (rmol * t_leaf(l))

END DO

! Transpiration can't be negative
el_sample(:,:) = MAX(0.0, el_sample(:,:))

! Approximate the leaf water potential and conductance.
CALL leaf_psi_jls( pft,                                                        &
                   som_n_sample,                                               &
                   land_pts,                                                   &
                   open_pts,                                                   &
                   veg_index,                                                  &
                   open_index,                                                 &
                   el_sample,                                                  &
                   psi_root_zone,                                              &
                 ! INTENT OUT
                   psi_sample,                                                 &
                   kl_sample                                                   &
  )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE stom_opt_mod_ci

END MODULE stom_opt_jls_mod