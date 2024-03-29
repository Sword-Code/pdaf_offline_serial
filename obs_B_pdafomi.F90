!$Id: obs_B_pdafomi.F90 793 2021-05-26 17:40:18Z lnerger $
!> PDAF-OMI observation module for type B observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!! TYPE = B
!!
!! __Observation type B:__
!! The observation type B in this tutorial are 6 observations at specified 
!! model grid points.
!!
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF.
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs)
!! and for the observation operator (obs_op).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!!
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! In addition, there are two optional routine, which are required if filters 
!! with localization are used:
!! * init_dim_obs_l_TYPE \n
!!           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays.
!! * localize_covar_TYPE \n
!!           Only required if the localized EnKF is used:
!!           Apply covariance localization in the LEnKF.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE obs_B_pdafomi

  USE mod_parallel, &
       ONLY: mype_filter    ! Rank of filter process
  USE PDAFomi, &
       ONLY: obs_f, obs_l   ! Declaration of observation data types
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_B        !< Whether to assimilate this data type
  !REAL    :: rms_obs_B      !< Observation error standard deviation (for constant errors)

  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.

    character(len=200), allocatable :: obsnames(:) !names of the observed variables
    integer :: nobsvar !number of argo variables
    integer, allocatable :: obsvar_p(:) !var number

! ***********************************************************************
! *** The following two data types are used in PDAFomi                ***
! *** They are declared in PDAFomi and only listed here for reference ***
! ***********************************************************************

! Data type to define the full observations by internally shared variables of the module
!   TYPE obs_f
!           Mandatory variables to be set in INIT_DIM_OBS
!      INTEGER :: doassim                   ! Whether to assimilate this observation type
!      INTEGER :: disttype                  ! Type of distance computation to use for localization
!                                           ! (0) Cartesian, (1) Cartesian periodic
!                                           ! (2) simplified geographic, (3) geographic haversine function
!      INTEGER :: ncoord                    ! Number of coordinates use for distance computation
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! Indices of observed field in state vector (process-local)
!           
!           Optional variables - they can be set in INIT_DIM_OBS
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!
!           Variables with predefined values - they can be changed in INIT_DIM_OBS
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!      INTEGER :: use_global_obs=1          ! Whether to use (1) global full obs. 
!                                           ! or (0) obs. restricted to those relevant for a process domain
!
!           The following variables are set in the routine PDAFomi_gather_obs
!      INTEGER :: dim_obs_p                 ! number of PE-local observations
!      INTEGER :: dim_obs_f                 ! number of full observations
!      INTEGER :: dim_obs_g                 ! global number of observations
!      INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
!      INTEGER :: off_obs_g                 ! Offset of this observation in overall global obs. vector
!      INTEGER :: obsid                     ! Index of observation over all assimilated observations
!      REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
!      REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
!      REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
!      INTEGER, ALLOCATABLE :: id_obs_f_lim(:) ! Indices of domain-relevant full obs. in global vector of obs.
!                                           ! (only if full obs. are restricted to process domain))
!   END TYPE obs_f

! Data type to define the local observations by internally shared variables of the module
!   TYPE obs_l
!      INTEGER :: dim_obs_l                 ! number of local observations
!      INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
!      INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector 
!      REAL, ALLOCATABLE :: distance_l(:)   ! Distances of local observations
!      REAL, ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
!      INTEGER :: locweight                 ! Specify localization function
!      REAL :: lradius                      ! localization radius
!      REAL :: sradius                      ! support radius for localization function
!   END TYPE obs_l
! ***********************************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could renamed the variables
  TYPE(obs_f), TARGET, PUBLIC :: thisobs      ! full observation
  TYPE(obs_l), TARGET, PUBLIC :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

CONTAINS

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations 
!! required for the analyses in the loop over all local 
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * thisobs\%doassim     - Whether to assimilate this type of observations
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!! * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
  SUBROUTINE init_dim_obs_B(step, dim_obs)

    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs
    USE mod_assimilation, &
         ONLY: nx, ny, nz, filtertype, local_range, varindex

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i, j, k, h                      ! Counters
    INTEGER :: cnt, cnt0                 ! Counters
    INTEGER :: dim_obs_p                 ! Number of process-local observations
    REAL, ALLOCATABLE :: obs_field(:,:,:,:)  ! Observation field read from file
    REAL, ALLOCATABLE :: std_field(:,:,:,:)  ! Observation std field read from file
    REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
    REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates 
    CHARACTER(len=2) :: stepstr          ! String for time step
    integer, parameter :: read_unit=12
    integer :: ios
    character(len=200) :: obsname !names of the observed variable

! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - obs type B'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_B) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2
    
    IF (ALLOCATED(obsnames)) DEALLOCATE(obsnames)
    
    open(unit=read_unit, file='data/obs/argo_names.txt', status='old', iostat=ios)
        if (ios /= 0) stop "Error opening file argo_names.txt"
        
        ios=0
        nobsvar=0
        do while (ios==0)
            read(read_unit, *, iostat=ios) obsname
            nobsvar=nobsvar+1
        end do
        nobsvar=nobsvar-1
        
        allocate(obsnames(nobsvar))
        
        rewind(read_unit)        
        read(read_unit, *, iostat=ios) obsnames
        
        do i=1,nobsvar
            write(*,*) trim(obsnames(i))
        end do

    close(read_unit)
    if (ios /= 0) stop "Error reading file argo_names.txt"

! **********************************
! *** Read PE-local observations ***
! **********************************

    ! Read observation field from file
    ALLOCATE(obs_field(nz, ny, nx, nobsvar))
    ALLOCATE(std_field(nz, ny, nx, nobsvar))

    IF (step<10) THEN
       WRITE (stepstr, '(i1)') step
    ELSE
       WRITE (stepstr, '(i2)') step
    END IF

    OPEN (read_unit, file='data/obs/argo.txt', status='old', iostat=ios)
        if (ios /= 0) stop "Error opening file argo.txt"
        READ (read_unit, *, iostat=ios) obs_field
    CLOSE (read_unit)
    if (ios /= 0) stop "Error reading file argo.txt"
    
    OPEN (read_unit, file='data/obs/argo_std.txt', status='old', iostat=ios)
        if (ios /= 0) stop "Error opening file argo_std.txt"
        READ (read_unit, *, iostat=ios) std_field
    CLOSE (read_unit)
    if (ios /= 0) stop "Error reading file argo_std.txt"


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***

    cnt = 0
    do h=1,nobsvar
        obsname=obsnames(h)
        if ((trim(obsname)=="Chl").or.(varindex(obsname)>0)) then
            DO j = 1, nx
                DO i= 1, ny
                    do k=1,nz
                        IF (obs_field(k,i,j,h) > -999.0) cnt = cnt + 1
                    end do
                END DO
            END DO
        end if
    end do
    dim_obs_p = cnt
    dim_obs = cnt

    IF (mype_filter==0) &
         WRITE (*,'(8x, a, i6)') '--- number of full observations', dim_obs
         
    IF (ALLOCATED(obsvar_p)) DEALLOCATE(obsvar_p)
    allocate(obsvar_p(dim_obs_p))


    ! *** Initialize vector of observations on the process sub-domain ***
    ! *** Initialize coordinate array of observations on the process sub-domain ***

    ! Allocate process-local observation arrays
    ALLOCATE(obs_p(dim_obs_p))
    ALLOCATE(ivar_obs_p(dim_obs_p))
    ALLOCATE(ocoord_p(3, dim_obs_p)) !probabilmente dovrei deifinire solo 2 righe e non 3, perche' ncoord=2

    ! Allocate process-local index array
    ! This array has a many rows as required for the observation operator
    ! 1 if observations are at grid points; >1 if interpolation is required
    ALLOCATE(thisobs%id_obs_p(1, dim_obs_p))

    cnt = 0
    do h=1,nobsvar
        obsname=obsnames(h)
        if ((trim(obsname)=="Chl").or.(varindex(obsname)>0)) then
            cnt0 = 0
            DO j = 1, nx
                DO i= 1, ny
                    do k=1, nz
                        cnt0 = cnt0 + 1
                        IF (obs_field(k, i, j, h) > -999.0) THEN
                            cnt = cnt + 1
                            thisobs%id_obs_p(1, cnt) = cnt0
                            obs_p(cnt) = obs_field(k, i, j, h)
                            ivar_obs_p(cnt) = 1.0 / std_field(k, i, j, h)**2
                            obsvar_p(cnt) = h
                            ocoord_p(1, cnt) = REAL(j)
                            ocoord_p(2, cnt) = REAL(i)
                            ocoord_p(3, cnt) = REAL(k) !probabilmente questo non serve, perche' ncoord=2, il che significa (credo) che la distanza viene calcolata solo sull'orizzontale. 
                        END IF
                    end do
                END DO
            END DO
        end if
    end do

! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ! *** Set inverse observation error variances ***

    !ivar_obs_p(:) = 1.0 / (rms_obs_B*rms_obs_B)


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, local_range, dim_obs)


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

!     IF (twin_experiment .AND. filtertype/=11) THEN
!        CALL read_syn_obs(file_syntobs_TYPE, dim_obs, thisobs%obs_f, 0, 1-mype_filter)
!     END IF


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    DEALLOCATE(obs_field)
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE init_dim_obs_B



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE obs_op_B(dim_p, dim_obs, state_p, ostate)
         
    USE mod_assimilation, &
        ONLY: nx, ny, nz, varindex
         
    USE PDAFomi_obs_f, & 
        ONLY: PDAFomi_gather_obsstate

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: ostate(dim_obs)       !< Full observed state

! *** Local variables ***
    
   integer :: i, idx
   real, ALLOCATABLE :: ostate_p(:)
   character(len=200) :: obsname !names of the observed variable
   

! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    IF (thisobs%doassim==1) THEN
       ! observation operator for observed grid point values
       
        ALLOCATE(ostate_p(thisobs%dim_obs_p))   
       
        do i=1, thisobs%dim_obs_p
        
            idx=thisobs%id_obs_p(1, i)
            obsname=obsnames(obsvar_p(i))
            
            if (trim(obsname)=="Chl") then
            
                ostate_p(i)=state_p(idx + (varindex("P1_Chl")-1)*nx*ny*nz) + &
                            state_p(idx + (varindex("P2_Chl")-1)*nx*ny*nz) + &
                            state_p(idx + (varindex("P3_Chl")-1)*nx*ny*nz) + &
                            state_p(idx + (varindex("P4_Chl")-1)*nx*ny*nz)
                            
            elseif (varindex(obsname)>0) then
            
                ostate_p(i)=state_p(idx + (varindex(obsname)-1)*nx*ny*nz)
                
            else 
                
                stop "something went wrong in obs_op_B"
                
            end if
            
        end do
        
        ! *** Global: Gather full observed state vector
        CALL PDAFomi_gather_obsstate(thisobs, ostate_p, ostate)

        ! *** Clean up
        DEALLOCATE(ostate_p)
       
    END IF

  END SUBROUTINE obs_op_B



!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
  SUBROUTINE init_dim_obs_l_B(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l

    ! Include localization radius and local coordinates
    USE mod_assimilation, &   
         ONLY: coords_l, local_range, locweight, srange

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step         !< Current time step
    INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
    INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
         locweight, local_range, srange, dim_obs_l)

  END SUBROUTINE init_dim_obs_l_B



!-------------------------------------------------------------------------------
!> Perform covariance localization for local EnKF on the module-type observation
!!
!! The routine is called in the analysis step of the localized
!! EnKF. It has to apply localization to the two matrices
!! HP and HPH of the analysis step for the module-type
!! observation.
!!
!! This routine calls the routine PDAFomi_localize_covar
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type.
!!
  SUBROUTINE localize_covar_B(dim_p, dim_obs, HP_p, HPH, coords_p)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_localize_covar

    ! Include localization radius and local coordinates
    USE mod_assimilation, &   
         ONLY: local_range, locweight, srange

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
    REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
    REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
    REAL, INTENT(in)    :: coords_p(:,:)         !< Coordinates of state vector elements


! *************************************
! *** Apply covariance localization ***
! *************************************

    CALL PDAFomi_localize_covar(thisobs, dim_p, locweight, local_range, srange, &
         coords_p, HP_p, HPH)

  END SUBROUTINE localize_covar_B



!-------------------------------------------------------------------------------
!> Implementation of adjoint observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE obs_op_adj_B(dim_p, dim_obs, ostate, state_p)

    USE mod_assimilation, &
         ONLY: nx, ny, nz, varindex

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: ostate(dim_obs)       !< Full observed state
    REAL, INTENT(inout) :: state_p(dim_p)        !< PE-local model state

! *** Local variables ***
    
    integer :: i, idx
    character(len=200) :: obsname !names of the observed variable

! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    IF (thisobs%doassim==1) THEN
       ! adjoint observation operator for observed grid point values
       
        state_p=0.0
        do i=1, thisobs%dim_obs_p 
            
            idx=thisobs%id_obs_p(1, i)
            obsname=obsnames(obsvar_p(i))
            
            if (trim(obsname)=="Chl") then
            
                state_p(idx + (varindex("P1_Chl")-1)*nx*ny*nz)=state_p(idx + (varindex("P1_Chl")-1)*nx*ny*nz)+ostate(thisobs%off_obs_f+i)
                state_p(idx + (varindex("P2_Chl")-1)*nx*ny*nz)=state_p(idx + (varindex("P2_Chl")-1)*nx*ny*nz)+ostate(thisobs%off_obs_f+i)
                state_p(idx + (varindex("P3_Chl")-1)*nx*ny*nz)=state_p(idx + (varindex("P3_Chl")-1)*nx*ny*nz)+ostate(thisobs%off_obs_f+i)
                state_p(idx + (varindex("P4_Chl")-1)*nx*ny*nz)=state_p(idx + (varindex("P4_Chl")-1)*nx*ny*nz)+ostate(thisobs%off_obs_f+i)
                            
            elseif (varindex(obsname)>0) then
            
                state_p(idx + (varindex(obsname)-1)*nx*ny*nz)=state_p(idx + (varindex(obsname)-1)*nx*ny*nz)+ostate(thisobs%off_obs_f+i)
                
            else
            
                stop "something went wrong in obs_op_adj_B"
            
            end if
            
        end do
       
    END IF

  END SUBROUTINE obs_op_adj_B

END MODULE obs_B_pdafomi
