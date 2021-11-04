!$Id: cvt_pdaf.F90 793 2021-05-26 17:40:18Z lnerger $
!> Apply covariance operator to a control vector
!!
!! The routine is called during the analysis step of
!! 3D-Var or hybrid 3D-Var. It has to apply the 
!! covariance operator (square root of P) to a vector 
!! in control space.
!!
!! For domain decomposition, the action is on
!! the control vector for the PE-local part of
!! the sub-state vector for the PE-local domain.
!!
!! This code variant uses an explicit array holding
!! the covariance operator as a matrix.
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE cvt_pdaf(iter, dim_p, dim_cvec, v_p, Vv_p)

  USE mod_assimilation, &
       ONLY: Vmat_p, nx, ny, nz, nvar, dim_eof_p, state3dvar_p

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: iter          !< Iteration of optimization
  INTEGER, INTENT(in) :: dim_p         !< PE-local observation dimension !!state dimension?
  INTEGER, INTENT(in) :: dim_cvec      !< Dimension of control vector
  REAL, INTENT(in)    :: v_p(dim_cvec) !< PE-local model state
  REAL, INTENT(inout) :: Vv_p(dim_p)   !< PE-local result vector
  
! **** internal variables ****
    real, allocatable :: Mv(:) !result of vertical operator
    real, allocatable :: totalchl(:,:,:)          !total chlorophyll
    real, allocatable :: field(:,:,:), field2(:,:,:,:) !chl field
    integer :: i, j, k
    
    allocate(Mv(dim_eof_p))
    allocate(field(nz,ny,nx))
    allocate(field2(nz,ny,nx,nvar))
    allocate(totalchl(nz,ny,nx))

! *********************
! *** Compute vertical operator ***
! *********************

  ! Transform control variable to state increment
  CALL dgemv('n', dim_eof_p, dim_cvec, 1.0, Vmat_p, &
       dim_eof_p, v_p, 1, 0.0, Mv, 1)
       
       

    DO j = 1, nx
        do i=1,ny
            field(1:nz,i, j) =   Mv(1 + (j-1)*ny*nz + (i-1)*nz : (j-1)*ny*nz + i*nz)
        end do
    END DO
    
! *********************
! *** Compute biogeochemical operator ***
! *********************

    totalchl=state3dvar_p(:,:,:,4)+state3dvar_p(:,:,:,9)+state3dvar_p(:,:,:,13)+state3dvar_p(:,:,:,17)
    
    do i=1,nvar
        field2(:,:,:,i)=field/totalchl*state3dvar_p(:,:,:,i)
    end do
    
    do k=1,nvar
        DO j = 1, nx
            do i=1, ny
                Vv_p(1 + (k-1)*nx*ny*nz + (j-1)*ny*nz + (i-1)*nz : & 
                    (k-1)*nx*ny*nz + (j-1)*ny*nz + i*nz) = field2(1:nz,i, j, k)            
            end do
        END DO
    end do
    
 ! *********************
! *** deallocate ***
! *********************   
    deallocate(Mv)
    deallocate(field)
    deallocate(field2)
   deallocate(totalchl)
   
   !write(*,*) "cvt v: ", v_p
   !write(*,*) "cvt Vv: ", Vv_p

END SUBROUTINE cvt_pdaf
