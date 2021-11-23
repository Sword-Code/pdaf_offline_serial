!$Id: initialize.F90 793 2021-05-26 17:40:18Z lnerger $
!BOP
!
! !ROUTINE: initialize  --- initialize the 2D offline example for PDAF
!
! !INTERFACE:
SUBROUTINE initialize()

! !DESCRIPTION:
! Routine to perform initialization of the 2D offline example for
! PDAF. Here, only the global size of the model domain and the
! global size of the model state vector need to be initialized.
! Generally, this could also be joined with the routine init_pdaf().
!
! For the 2D offline tutorial example, the domain is defined by the
! dimensions nx and ny. The state vector size is nx*ny.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, & ! Model variables
       ONLY: dim_state_p, nx, ny, nz, nvar, dim_eof_p, varnames
  USE mod_parallel, &     ! Parallelization variables
       ONLY: MPI_COMM_WORLD, init_parallel, finalize_parallel

    IMPLICIT NONE
    integer :: i , ios
    integer, parameter :: read_unit = 11

!EOP

! *** Model specifications ***
    nx = 1 !36    ! Extent of grid in x-direction
    ny = 1 !18    ! Extent of grid in y-direction
    nz = 196 !10    ! Extent of grid in z-direction
    nvar = 17 !2   ! number of variables
    
    dim_state_p   = nx * ny * nz * nvar ! State dimension (shared via MOD_OFFLINE)
    dim_eof_p = nx * ny * nz
    
    allocate(varnames(nvar))
    
    open(unit=read_unit, file='data/init/names.txt', status='old', iostat=ios)

        do i = 1, nvar
            read(read_unit, *, iostat=ios) varnames(i)
            if (ios /= 0) stop "Error reading file names.txt"
        end do

    close(read_unit)


! *** Screen output ***
    WRITE (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
    WRITE (*, '(22x,a)') 'MODEL: 2D Offline Example for Tutorial'
    WRITE (*, *) 'Grid size:',nx,'x',ny,'x',nz
    WRITE (*, *) 'Number of variables:',nvar
    WRITE (*, *) 'Variable list:'
  
    do i = 1, nvar
        WRITE (*, *) trim(varnames(i))
    end do
    
    WRITE (*, '(5x, a, i7)') &
        'Global model state dimension:', dim_state_p

END SUBROUTINE initialize
