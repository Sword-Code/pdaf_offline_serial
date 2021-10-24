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
       ONLY: dim_state_p, nx, ny, nz
  USE mod_parallel, &     ! Parallelization variables
       ONLY: MPI_COMM_WORLD, init_parallel, finalize_parallel

  IMPLICIT NONE

!EOP

! *** Model specifications ***
  nx = 36    ! Extent of grid in x-direction
  ny = 18    ! Extent of grid in y-direction
  nz = 10    ! Extent of grid in z-direction

  dim_state_p   = nx * ny * nz ! State dimension (shared via MOD_OFFLINE)



! *** Screen output ***
  WRITE (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
  WRITE (*, '(22x,a)') 'MODEL: 2D Offline Example for Tutorial'
  WRITE (*, *) 'Grid size:',nx,'x',ny,'x',nz
  WRITE (*, '(5x, a, i7)') &
       'Global model state dimension:', dim_state_p

END SUBROUTINE initialize
