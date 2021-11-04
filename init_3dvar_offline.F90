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
       ONLY: nx, ny, nz, nvar, dim_cvec, Vmat_p, dim_eof_p, state3dvar_p

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
  INTEGER :: i, j, k, member  ! Counters
  REAL, ALLOCATABLE :: field(:,:,:,:)     ! global model field
  REAL, ALLOCATABLE :: eof(:)     ! eof
  REAL, ALLOCATABLE :: eoffield(:,:,:)     ! eof field

! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Read initial state and generate square root of B ***
  ! *** by reading the full ensemble on filter-PE 0      ***
  WRITE (*, '(/9x, a)') 'Initialize state and B^1/2 for 3D-Var'
  WRITE (*, '(9x, a)') '--- read ensemble from files'
  WRITE (*, '(9x, a, i5)') '--- members in B^1/2:  ', dim_cvec
  

  ! allocate memory for temporary fields
  ALLOCATE(field(nz, ny, nx, nvar))
  ALLOCATE(eof(nz))
  ALLOCATE(eoffield(nz, ny, nx))

  ! Allocate matrix holding B^1/2 (from mod_assimilation)
  ALLOCATE(Vmat_p(dim_eof_p, dim_cvec))

  ! Allocate initial state
  ALLOCATE(state3dvar_p(nz, ny, nx, nvar))


! **********************************************
! *** Initialize square-root of P for 3D-Var ***
! **********************************************

  WRITE (*, '(9x, a)') 'Initialize B^1/2'
  
    OPEN(11, file = 'data/eof.txt', status='old')
        DO member = 1, dim_cvec
            READ (11, *) eof(:)
            
            do i=1,nx
                do j=1,ny
                    eoffield(:,j,i)=eof
                end do
            end do
            
            do j=1,nx
                do i=1,ny
                    Vmat_p(1 + (j-1)*ny*nz + (i-1)*nz : & 
                        (j-1)*ny*nz + i*nz, member) = eoffield(1:nz,i, j)            
                end do
            end do
        end do
    
    CLOSE(11)


! ******************************************
! *** Initialize ensemble array for PDAF ***
! ******************************************

    OPEN(11, file = 'data/phyto.txt', status='old')
        DO member = 1, nvar
            READ (11, *) eof(:)
            
            do i=1,nx
                do j=1,ny
                    field(:,j,i,member)=eof
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
    
    ens_p(:,1) = state_p(:)
    state3dvar_p=field
    
    OPEN(11, file = 'chl_init.txt', status = 'replace')
 
        WRITE (11, *) field(:,1, 1,4) +field(:,1, 1,9) +field(:,1, 1, 13) +field(:,1, 1,17)

    CLOSE(11)


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(field)
  DEALLOCATE(eof)
  DEALLOCATE(eoffield)

END SUBROUTINE init_3dvar_offline
