!$Id: cvt_adj_pdaf.F90 793 2021-05-26 17:40:18Z lnerger $
!> Apply adjoint covariance operator to a state vector
!!
!! The routine is called during the analysis step.
!! It has to apply the adjoint covariance operator 
!! (transpose of square root of P) to a vector in
!! state space.
!!
!! For domain decomposition, the action is for
!! the PE-local sub-domain of the state. Thus the
!! covariance operator is applied to a sub-state.
!!
!! This code variant uses an explicit array holding
!! the covariance operator as a matrix.
!!
!! __Revision history:__
!! * 2021-03 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE cvt_adj_pdaf(iter, dim_p, dim_cvec, Vv_p, v_p)

  USE mod_assimilation, &
       ONLY: Vmat_p, nx, ny, nz, nvar, dim_eof_p, state3dvar_p, varindex

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: iter          !< Iteration of optimization
  INTEGER, INTENT(in) :: dim_p         !< PE-local observation dimension
  INTEGER, INTENT(in) :: dim_cvec      !< Dimension of control vector
  REAL, INTENT(in)    :: Vv_p(dim_p)   !< PE-local input vector
  REAL, INTENT(inout) :: v_p(dim_cvec) !< PE-local result vector
  
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
! *** Compute biogeochemical operator ***
! *********************

    totalchl=state3dvar_p(:,:,:,varindex("P1_Chl"))+state3dvar_p(:,:,:,varindex("P2_Chl"))+ &
                state3dvar_p(:,:,:,varindex("P3_Chl"))+state3dvar_p(:,:,:,varindex("P4_Chl")) 
    
    do k=1, nvar
        DO j = 1, nx
            do i=1,ny
                field2(1:nz,i, j, k) =   Vv_p(1 + (k-1)*nx*ny*nz + (j-1)*ny*nz + (i-1)*nz : & 
                                        (k-1)*nx*ny*nz + (j-1)*ny*nz + i*nz)
            end do
        END DO
     end do
    
    field=0.0
    do i=1,nvar
        field=field + field2(:,:,:,i)*state3dvar_p(:,:,:,i)
    end do
    
    field=field/totalchl
    
    
    DO j = 1, nx
        do i=1, ny
            Mv(1 + (j-1)*ny*nz + (i-1)*nz : (j-1)*ny*nz + i*nz) = field(1:nz,i, j)            
        end do
    END DO
    


! ***********************
! *** Compute vertical operator ***
! ***********************

  ! Transform control variable to state increment
  CALL dgemv('t', dim_eof_p, dim_cvec, 1.0, Vmat_p, &
       dim_eof_p, Mv, 1, 0.0, v_p, 1)
       
    
 ! *********************
! *** deallocate ***
! *********************   
    deallocate(Mv)
    deallocate(field)
    deallocate(field2)
   deallocate(totalchl)
   
   !write(*,*) "adj Vv: ",  Vv_p
   !write(*,*) "adj v: ",  v_p

END SUBROUTINE cvt_adj_pdaf
