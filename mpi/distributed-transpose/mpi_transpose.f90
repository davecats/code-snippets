!============================================!
!                                            !
!             Distributed Traspose           !
!                  for the                   !
!      Direct Numerical Simulation (DNS)     !
!        of a turbulent channel flow         !
!                                            !
!============================================!
! 
! Author: M.Sc. Davide Gatti
! Date  : 9/Sep/2015
! 

MODULE mpi_transpose

  USE, intrinsic :: iso_c_binding
  USE typedef
  IMPLICIT NONE

CONTAINS

  ! First transpose
  ! (nx0:nxN)-(0:nzd) --> (nz0:nzN)-(0:nx)
  ! Locally transform and dot product
  ! (nz0:nzN)-(0:nx)  --> (nz0:nzN)-(0:nxd)
  ! Re-transform back and copy
  ! (nz0:nzN)-(0:nx)  <-- (nz0:nzN)-(0:nxd)
  ! Transpose back
  ! (nx0:nxN)-(0:nzd) <-- (nz0:nzN)-(0:nx)

  ! We need space for max of (nx0:nxN)-(0:nzd) --> (nz0:nzN)-(0:nx)
  ! and for the transform (nz0:nzN)-(0:nxd)

  !--- Find size to allocate to calling process ---! 
  !------------------------------------------------!
  SUBROUTINE problem_size(iproc,nproc,nx,nz,nx0,nxN,nxB,nz0,nzN,nzB,block)
    integer(C_INT), intent(in)  :: nx,nz,iproc,nproc
    integer(C_INT), intent(out) :: nx0,nxN,nxB,nz0,nzN,nzB,block
    nx0=iproc*nx/nproc;   nxN=(iproc+1)*nx/nproc-1;  nxB=nxN-nx0+1;
    nz0=iproc*nz/nproc;   nzN=(iproc+1)*nz/nproc-1;  nzB=nzN-nz0+1;
    block=max(nxB*nz,nx*nzB)
    WRITE(*,*) "iproc=",iproc," nx0=",nx0," nxN=",nxN," nxB=",nxB, "nz0=",nz0," nzN=",nzN," nzB=",nzB, " s=", block 
  END SUBROUTINE problem_size
  

END MODULE mpi_transpose
