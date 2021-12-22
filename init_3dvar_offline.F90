!$Id: init_3dvar_offline.F90 825 2021-10-13 10:41:48Z lnerger $
!>  Initialize 3D-Var
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in 3D-Var
!!
!! The routine is called when the filter is
!! initialized in PDAF_filter_init.  It has
!! to initialize an ensemble of dim_ens states.
!!
!! The routine is called by all filter processes and 
!! initializes the ensemble for the PE-local domain.
!!
!! Implementation for the 2D online example
!! without parallelization. Here, the ensemble 
!! information is directly read from files.
!!
!! __Revision history:__
!! * 2021-05 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_3dvar_offline(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

  USE mod_assimilation, &
       ONLY: nx, ny, nz, nvar, state3dvar_p, varindex

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: filtertype                !< Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                     !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                   !< Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)            !< PE-local model state
  !< (It is not necessary to initialize the array 'state_p' for ensemble filters.
  !< It is available here only for convenience and can be used freely.)
  REAL, INTENT(inout) :: Uinv(1,1)                 !< Array not referenced for 3D-Var
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)     !< PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                   !< PDAF status flag

! *** local variables ***
  INTEGER :: i, j, k, member, ios  ! Counters
  REAL, ALLOCATABLE :: field(:,:,:,:)     ! global model field
  REAL, ALLOCATABLE :: column(:)     ! tracer in water column
  REAL :: invdim_ens                   ! Inverse ensemble size
  CHARACTER(len=2) :: ensstr          ! String for ensemble member

! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Read initial state and generate square root of B ***
  ! *** by reading the full ensemble on filter-PE 0      ***
  WRITE (*, '(/9x, a)') 'Initialize state and B^1/2 for 3D-Var'
  

  ! allocate memory for temporary fields
  ALLOCATE(field(nz, ny, nx, nvar))
  ALLOCATE(column(nz))

  ! Allocate initial state
  ALLOCATE(state3dvar_p(nz, ny, nx, nvar))


! **********************************************
! *** Initialize square-root of P for 3D-Var ***
! **********************************************

    call read_eof_3dvar


! ******************************************
! *** Initialize ensemble array for PDAF ***
! ******************************************

    WRITE (*, '(/9x, a)') 'Initialize state for 3D-Var'
    
    if (.false.) then

        OPEN(11, file = 'data/forecast/phyto.txt', status='old')
            DO member = 1, nvar
                READ (11, *) column(:)
                
                do i=1,nx
                    do j=1,ny
                        field(:,j,i,member)=column
                    end do
                end do
            end do
        CLOSE(11)
                
        do k=1,nvar
            DO j = 1, nx
                do i=1, ny
                    state_p(1 + (k-1)*nx*ny*nz + (j-1)*ny*nz + (i-1)*nz : & 
                        (k-1)*nx*ny*nz + (j-1)*ny*nz + i*nz) = field(1:nz,i, j, k)            
                end do
            END DO
        end do
        
    else
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
        
        field=reshape(state_p, (/nz, ny, nx, nvar/))
        
    end if
    
    ens_p(:,1) = state_p(:)
    state3dvar_p=field
    
    OPEN(11, file = 'data/diag/chl_init.txt', status = 'replace')
 
        WRITE (11, *) field(:,1, 1,varindex("P1_Chl")) +field(:,1, 1,varindex("P2_Chl")) +field(:,1, 1, varindex("P3_Chl")) +field(:,1, 1,varindex("P4_Chl"))

    CLOSE(11)


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(field)
  DEALLOCATE(column)


END SUBROUTINE init_3dvar_offline
