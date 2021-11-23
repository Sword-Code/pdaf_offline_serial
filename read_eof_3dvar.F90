 
!! code by Simone Spada
SUBROUTINE read_eof_3dvar

  USE mod_assimilation, &
       ONLY: nx, ny, nz, dim_cvec, Vmat_p, dim_eof_p

  IMPLICIT NONE

! *** local variables ***
  INTEGER :: i, j, member  ! Counters
  REAL, ALLOCATABLE :: eof(:)     ! eof
  REAL, ALLOCATABLE :: eoffield(:,:,:)     ! eof field

! **********************
! *** INITIALIZATION ***
! **********************

  ! *** Read initial state and generate square root of B ***
  ! *** by reading the full ensemble on filter-PE 0      ***
  WRITE (*, '(/9x, a)') 'Initialize B^1/2 for 3D-Var'
  WRITE (*, '(9x, a)') '--- read eofs from files'
  WRITE (*, '(9x, a, i5)') '--- eofs in B^1/2:  ', dim_cvec
  

  ! allocate memory for temporary fields
  ALLOCATE(eof(nz))
  ALLOCATE(eoffield(nz, ny, nx))

  ! Allocate matrix holding B^1/2 (from mod_assimilation)
  ALLOCATE(Vmat_p(dim_eof_p, dim_cvec))


! **********************************************
! *** Initialize square-root of P for 3D-Var ***
! **********************************************

    OPEN(11, file = 'data/init/eof.txt', status='old')
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


! ****************
! *** clean up ***
! ****************

  DEALLOCATE(eof)
  DEALLOCATE(eoffield)

END SUBROUTINE read_eof_3dvar
