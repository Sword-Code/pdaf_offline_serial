!$Id: init_ens_offline.F90 793 2021-05-26 17:40:18Z lnerger $
!BOP
!
! !ROUTINE: init_ens_offline --- Initialize ensemble for SEIK in offline mode
!
! !INTERFACE:
SUBROUTINE init_ens_offline(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize an ensemble of dim\_ens states.
! For the offline mode, the ensemble will be
! typically read-in from files.
!
! The routine is called by all filter processes and 
! initializes the ensemble for the PE-local domain.
!
! Implementation for the 2D offline example
! without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: nx, ny, nz, nvar, state3dvar_p, varindex

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  ! It is not necessary to initialize the array 'state_p' for SEIK. 
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init    (as U_ens_init)
!EOP

! *** local variables ***
  INTEGER :: member, ios             ! Counters
  CHARACTER(len=2) :: ensstr          ! String for ensemble member
  REAL :: invdim_ens                   ! Inverse ensemble size


! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Generate full ensemble on filter-PE 0 ***
  WRITE (*, '(/9x, a)') 'Initialize state ensemble'
  WRITE (*, '(9x, a)') '--- read ensemble from files'
  WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens
  

  ! Allocate initial state
  ALLOCATE(state3dvar_p(nz, ny, nx, nvar))


! ********************************
! *** Read ensemble from files ***
! ********************************
    
    state_p=0.0
    DO member = 1, dim_ens
        WRITE (ensstr, '(i2.2)') member
        OPEN(11, file = 'data/forecast/ens_'//TRIM(ensstr)//'.txt', status='old', iostat=ios)
            
            if (ios/=0) stop "Error opening forecast ensemble file"
            
            read(11,*) ens_p(:,member)

        CLOSE(11)
        
        state_p=state_p+ens_p(:,member)
    END DO
    
    invdim_ens    = 1.0 / REAL(dim_ens)
    state_p=state_p*invdim_ens
    
    state3dvar_p=reshape(state_p, (/nz, ny, nx, nvar/))
        
    OPEN(11, file = 'data/diag/chl_init.txt', status = 'replace')
 
        WRITE (11, *) state3dvar_p(:,1, 1,varindex("P1_Chl")) + state3dvar_p(:,1, 1,varindex("P2_Chl")) + & 
            state3dvar_p(:,1, 1, varindex("P3_Chl")) + state3dvar_p(:,1, 1,varindex("P4_Chl"))

    CLOSE(11)
    

! *****************************************************
! *** Initialize square-root of P for hybrid 3D-Var ***
! *****************************************************

  IF (filtertype==13) THEN
     
    call read_eof_3dvar

  END IF

END SUBROUTINE init_ens_offline
