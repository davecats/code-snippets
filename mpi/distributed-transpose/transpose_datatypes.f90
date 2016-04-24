!============================================!
!                                            !
!     Direct Numerical Simulation (DNS)      !
!       of a turbulent channel flow          !
!                                            !
!============================================!
! 
! Author: M.Sc. Davide Gatti
! Date  : 28/Jul/2015
! 

PROGRAM channel

  
  USE, intrinsic :: iso_c_binding
  USE mpi_f08
  USE mpi_transpose

  implicit none
  
  integer(C_INT) :: nproc,iproc,ierr
  integer(C_INT), parameter :: nx=8, nz=4, nxd=3*nx/2, nzd=3*nz
  integer(C_INT) :: nx0,nxN,nxB,nz0,nzN,nzB
  integer(C_INT) :: ix,iy,iz,block
  complex(C_DOUBLE_COMPLEX), allocatable :: AA(:,:,:),BB(:,:,:)
  TYPE(MPI_Datatype) :: row, matrix, column, trans, tmp, vel, cmpl
  TYPE(MPI_Status) :: stat
  integer(kind=MPI_ADDRESS_KIND) :: stride,lb

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  IF (iproc==0) THEN 
    WRITE(*,*) " "
    WRITE(*,*) "!====================================================!"
    WRITE(*,*) "!                 T R A N S P O S E                  !"
    WRITE(*,*) "!====================================================!"
  END IF

  CALL problem_size(iproc,nproc,nx,nzd,nx0,nxN,nxB,nz0,nzN,nzB,block)
  ALLOCATE(BB(1:6,0:nxd-1,nz0:nzN))
  ALLOCATE(AA(1:6,0:nzd-1,nx0:nxN))


  AA=0;

  CALL MPI_Type_contiguous(2,MPI_DOUBLE_PRECISION,cmpl); 
  CALL MPI_Type_contiguous(3,cmpl,vel); tmp=vel;
  lb=0; stride=6*16; CALL MPI_Type_create_resized(tmp,lb,stride,vel);

  CALL MPI_Type_vector(nxB,nzB,nzd,vel,row);
  lb=0; stride=6*16*nzB; CALL MPI_Type_create_resized(row,lb,stride,matrix);
  CALL MPI_Type_commit(row)    
  CALL MPI_Type_commit(matrix)    

  CALL MPI_Type_vector(nzB,1,nxd,vel,column)
  lb=0; stride=6*16;  CALL MPI_Type_create_resized(column,lb,stride,tmp)
  CALL MPI_Type_contiguous(nxB,tmp,trans);
  CALL MPI_Type_commit(column)   
  CALL MPI_Type_commit(trans)    
                           

  DO ix=nx0,nxN
    DO iz=0,nzd-1
      AA(1,iz,ix)=ix*(nzd)+iz
      AA(2,iz,ix)=2000+ix*(nzd)+iz
      AA(3,iz,ix)=3000+ix*(nzd)+iz
      AA(4,iz,ix)=4000+ix*(nzd)+iz
      AA(5,iz,ix)=5000+ix*(nzd)+iz
      AA(6,iz,ix)=6000+ix*(nzd)+iz
    END DO
  END DO 

  CALL MPI_Barrier(MPI_COMM_WORLD)

  DO ix=nx0,nxN
    WRITE(*,*) ix, int(AA(1,:,ix))
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD) 

  !Transpose externally
  CALL MPI_Alltoall(AA, 1, matrix, BB, 1, trans, MPI_COMM_WORLD)
  
  DO iz=nz0,nzN
    WRITE(*,*) iz, int(BB(1,:,iz))
  END DO 


  CALL MPI_Type_free(matrix);
  CALL MPI_Type_free(row);
  CALL MPI_Type_free(column);
  CALL MPI_Type_free(trans);
  CALL MPI_FINALIZE(ierr)

END PROGRAM channel
