!$Id: prepoststep_3dvar_offline.F90 793 2021-05-26 17:40:18Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_3dvar_offline --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_3dvar_offline(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
! 
! The routine is called for global filters (e.g. SEIK)
! before the analysis and after the ensemble transformation.
! For local filters (e.g. LSEIK) the routine is called
! before and after the loop over all local analysis
! domains.
! The routine provides full access to the state 
! estimate and the state ensemble to the user.
! Thus, user-controlled pre- and poststep 
! operations can be performed here. For example 
! the forecast and the analysis states and ensemble
! covariance matrix can be analyzed, e.g. by 
! computing the estimated variances. 
! For the offline mode, this routine is the place
! in which the writing of the analysis ensemble
! can be performed.
!
! If a user considers to perform adjustments to the 
! estimates (e.g. for balances), this routine is 
! the right place for it.
!
! Implementation for the 2D offline example
! without parallelization.
!
! !REVISION HISTORY:
! 2021-05 - Lars Nerger - Initial code based on prepoststep_ens_offline
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: nx, ny, nz, nvar, dim_cvec, Vmat_p, dim_eof_p, varindex

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step (not relevant for offline mode)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_seik_update    (as U_prepoststep)
! Called by: PDAF_lseik_update    (as U_prepoststep)
! Calls: PDAF_add_increment
! Calls: PDAF_seik_TtimesA
! Calls: memcount
! Calls: dgemm (BLAS)
! Calls: dgesv (LAPACK)
! Calls: MPI_send
! Calls: MPI_recv
!EOP

! *** local variables ***
  INTEGER :: i, j, k, member             ! counters
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  REAL :: invdim_ens                  ! Inverse ensemble size
!  REAL :: invdim_ensm1                ! Inverse of ensemble size minus 1
  REAL :: rmserror_est                ! estimated RMS error
  REAL, ALLOCATABLE :: variance(:)    ! model state variances
  REAL, ALLOCATABLE :: field(:,:,:,:)     ! global model field
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  REAL :: fact                        ! Scaling factor


! **********************
! *** INITIALIZATION ***
! **********************

  IF (firsttime) THEN
     WRITE (*, '(8x, a)') 'Analyze forecasted state for 3D-Var'
  ELSE
     WRITE (*, '(8x, a)') 'Analyze and write assimilated state for 3D-Var'
  END IF

  ! Allocate fields
  ALLOCATE(variance(dim_eof_p))

  ! Initialize numbers
  rmserror_est  = 0.0
  invdim_ens    = 1.0 / REAL(dim_ens)  
!  invdim_ensm1  = 1.0 / REAL(dim_ens - 1)

!*****************************************************************************************
!*** Cut off negative values ***
!*******************************************

    do i=1,dim_p
        if (ens_p(i,1)<0.0) then 
            ens_p(i,1)=0.0
        end if
    end do


! **************************************************************
! *** Perform prepoststep for 3D-Var in which dim_ens=1      ***
! *** The sampled error is here computed from B^(1/2)        ***
! **************************************************************

  ! *** Initialize state estimate  
  state_p(:) = ens_p(:,1)

  ! *** Compute sampled variances ***
  variance(:) = 0.0
  DO member = 1, dim_cvec
     DO j = 1, dim_eof_p
        variance(j) = variance(j) &
             + Vmat_p(j,member) * Vmat_p(j,member)
     END DO
  END DO


! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  ! total estimated RMS error
  DO i = 1, dim_eof_p
     rmserror_est = rmserror_est + variance(i)
  ENDDO
  rmserror_est = SQRT(rmserror_est / dim_eof_p)


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
  WRITE (*, '(12x, a, es12.4)') &
       'RMS error according to sampled variance: ', rmserror_est

  
! *******************
! *** File output ***
! *******************

  IF (.not. firsttime) THEN

     WRITE (*, '(8x, a)') '--- write ensemble and state estimate'

     ALLOCATE(field(nz, ny, nx, nvar))

     
if (.false.) then
     ! Write analysis ensemble
     DO member = 1, dim_ens
        do k=1,nvar
            DO j = 1, nx
                do i=1,ny
                    field(1:nz,i, j, k) =   ens_p(1 + (k-1)*nx*ny*nz + (j-1)*ny*nz + (i-1)*nz : & 
                                            (k-1)*nx*ny*nz + (j-1)*ny*nz + i*nz, member)
                end do
            END DO
        end do

        WRITE (ensstr, '(i2.2)') member
        OPEN(11, file = 'data/analysis/ens_'//TRIM(ensstr)//'_ana.txt', status = 'replace')
 
        DO i = 1, ny
           WRITE (11, *) field(1,i, :,1)
        END DO

        CLOSE(11)
     END DO
end if

     ! Write analysis state
     do k=1, nvar
        DO j = 1, nx
            do i=1,ny
                field(1:nz,i, j, k) =   state_p(1 + (k-1)*nx*ny*nz + (j-1)*ny*nz + (i-1)*nz : & 
                                        (k-1)*nx*ny*nz + (j-1)*ny*nz + i*nz)
            end do
        END DO
     end do
    
     OPEN(11, file = 'data/analysis/state_ana.txt', status = 'replace')
 
        DO i = 1, nvar
            WRITE (11, *) field(:,1, 1,i)
        END DO

     CLOSE(11)
     
     OPEN(11, file = 'data/diag/chl.txt', status = 'replace')
        
        WRITE (11, *) field(:,1, 1,varindex("P1_Chl")) +field(:,1, 1,varindex("P2_Chl")) +field(:,1, 1, varindex("P3_Chl")) +field(:,1, 1,varindex("P4_Chl"))

     CLOSE(11)


     DEALLOCATE(field)
     
  END IF


! ********************
! *** finishing up ***
! ********************

  DEALLOCATE(variance)

  firsttime = .FALSE.

END SUBROUTINE prepoststep_3dvar_offline
