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
 
  USE mpi_f08
  USE mpi_transpose
  
  integer(C_INT) :: nproc,iproc,ierr
  integer(C_INT) :: nx=8,nzd=4,nx0,nxN,nxB,nz0,nzN,nzB
  integer(C_INT) :: block,ix,iy,iz,i,j
  real(C_DOUBLE), allocatable, target :: AA(:,:),AAt(:,:),Ain(:),Aout(:)
  !real(C_DOUBLE), pointer :: AAt(:,:)

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
  ALLOCATE(AA(nx0:nxN,0:nzd-1),Ain(0:nxB*nzd-1))
  ALLOCATE(AAt(nz0:nzN,0:nx-1),Aout(0:nzB*nx-1))
  !AAt(0:nx-1,nz0:nzN) => AA

  CALL MPI_Barrier(MPI_COMM_WORLD)

  DO ix=nx0,nxN
    DO iz=0,nzd-1
      AA(ix,iz)=ix*(nzd)+iz
    END DO
  END DO 

  CALL MPI_Barrier(MPI_COMM_WORLD)

  DO ix=nx0,nxN
    WRITE(*,*) ix, int(AA(ix,:))
  END DO

   CALL MPI_Barrier(MPI_COMM_WORLD)
   WRITE(*,*) "  "
   WRITE(*,*) "  "
   WRITE(*,*) "  "
    CALL MPI_Barrier(MPI_COMM_WORLD)

  !Pack
  i=0 
  DO j=0,nproc-1
    DO ix=nx0,nxN
      Ain(i:i+nzd/nproc-1)=AA(ix,j*(nzd/nproc):(j+1)*(nzd/nproc)-1)
      i=i+nzd/nproc
    END DO
  END DO

  !Send
  CALL MPI_Alltoall(Ain, nxB*nzB, MPI_DOUBLE_PRECISION, Aout, nxB*nzB, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD)

  !Unpack
  i=0
  DO ix=0,nx-1
    DO iz=nz0,nzN
      AAt(iz,ix)=Aout(i)
      i=i+1
    END DO
  END DO

  DO iz=nz0,nzN
    WRITE(*,*) iz, int(AAt(iz,:))
  END DO

  CALL MPI_FINALIZE(ierr)

END PROGRAM channel
