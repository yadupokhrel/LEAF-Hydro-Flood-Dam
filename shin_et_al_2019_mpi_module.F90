!
!  Reservoir scheme of LEAF-Hydro-Flood-Dam (LHFD) model
!
!  Shin, S., Y. Pokhrel, and G. Miguez-Macho (2019) 
!     “High Resolution Modeling of Reservoir Release and Storage Dynamics at the Continental Scale”, 
!     Water Resources Research, Vol. 55, Issue 1, pp.787-810. doi:10.1029/2018WR023025
!
!  Originally authored by Gonzalo Miguez-Macho (https://scholar.google.com/citations?user=wSKBQLsAAAAJ).
!  Modified by Sanghoon Shin (https://sites.google.com/view/sshin) & Yadu Pokhrel (https://water.egr.msu.edu/)
!
!
MODULE module_parallel
    
	! Instead of "use mpi", "include 'mpif.h'" can also be used, e.g., when you use "Intel/15.0 & OpenMPI/1.10.0"
    use mpi
    use module_timing

    implicit none
     
    !include 'mpif.h' 

    integer, parameter :: npmax=300

    integer, save, dimension (0:npmax) :: domblock,domblock2dint,domblock3d,domblock3dint &
                            ,domblock3dnzs,domblocksmall,domblock3dsmall &
                            ,domblock3dnzssmall,nini,nend,rcountblock,rcountblocksmall,disp &
                            ,domblock312small
    integer, save :: columntype,columntype2
    integer, save :: pid,numtasks
    integer :: status(MPI_STATUS_SIZE),ierr

    integer, parameter :: n2big=1450,n3big=1510,nw=1,ne=n2big,ns=1,nn=n3big

! npmax		the maximumn number of processors
! n2big		total number of grid cells in x-direction
! n3big		total number of grid cells in y-direction

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INITIALIZEDOMAIN(n2,n3,nzg,nzs)
    integer :: n2,n3,nzg,nzs
    integer :: n,nmax
    integer :: tasktype
    integer :: request
    integer, allocatable :: req(:),stats(:,:)
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

    if(numtasks.eq.1)then
        nini(0)=1
        nend(0)=n3
        return
    endif

    rcountblocksmall(0)=n2
    rcountblock(0)=0
    disp(0)=0

    rcountblocksmall(numtasks-1)=n2
    rcountblock(numtasks-1)=0
    disp(numtasks-1)=0

    call MPI_TYPE_CONTIGUOUS(numtasks,MPI_INTEGER,tasktype,ierr)
    call MPI_Type_commit(tasktype,ierr)

    if(pid.eq.0)then

        nini(0)=0
        nini(numtasks-1)=n3+1
        call dividedomain(n2,n3,nini)

        allocate(req(numtasks-1))
        allocate(stats(MPI_STATUS_SIZE,numtasks-1))

        do n=1,numtasks-1
            call MPI_isend(nini(0),1,tasktype,n,1,MPI_COMM_WORLD,req(n),ierr)
        enddo
        if(numtasks.gt.1)call MPI_waitall(numtasks-1,req,stats,ierr)

        deallocate(req,stats)

    else
        call MPI_irecv(nini(0),1,tasktype,0,1,MPI_COMM_WORLD,request,ierr)
        call MPI_wait(request,status,ierr)
    endif

    call MPI_TYPE_FREE (tasktype,ierr)

    nend(0)=0
    nend(numtasks-1)=n3+1
    nend(numtasks-2)=n3

    do n=2,numtasks-2
        nend(n-1)=nini(n)+1
    enddo

    !gmmdeclare pieces to be send and received
    do n=1,numtasks-2

        nmax=nend(n)-nini(n)+1

        if(pid.eq.0)write(6,*)nini(n),nend(n),nmax,n,pid

        rcountblocksmall(n)=n2*(nmax-2)
        if(n.eq.1.or.n.eq.numtasks-1)rcountblocksmall(n)=n2*(nmax-1)
        rcountblock(n)=n2*nmax

        disp(n)=rcountblocksmall(n-1)+disp(n-1)

        call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblock(n),ierr)
        call MPI_Type_commit(domblock(n),ierr)
        call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_INTEGER,domblock2dint(n),ierr)
        call MPI_Type_commit(domblock2dint(n),ierr)
        call MPI_TYPE_CONTIGUOUS(nzg*n2*nmax,MPI_REAL,domblock3d(n),ierr)
        call MPI_Type_commit(domblock3d(n),ierr)
        call MPI_TYPE_CONTIGUOUS(nzs*n2*nmax,MPI_REAL,domblock3dnzs(n),ierr)
        call MPI_Type_commit(domblock3dnzs(n),ierr)
        call MPI_TYPE_CONTIGUOUS(nzg*n2*nmax,MPI_INTEGER,domblock3dint(n),ierr)
        call MPI_Type_commit(domblock3dint(n),ierr)
        call MPI_TYPE_CONTIGUOUS(rcountblocksmall(n),MPI_REAL,domblocksmall(n),ierr)
        call MPI_Type_commit(domblocksmall(n),ierr)
        call MPI_TYPE_CONTIGUOUS(nzg*rcountblocksmall(n),MPI_REAL,domblock3dsmall(n),ierr)
        call MPI_Type_commit(domblock3dsmall(n),ierr)
        call MPI_TYPE_CONTIGUOUS(nzs*rcountblocksmall(n),MPI_REAL,domblock3dnzssmall(n),ierr)
        call MPI_Type_commit(domblock3dnzssmall(n),ierr)

        call MPI_TYPE_CONTIGUOUS(12*n2*nmax,MPI_REAL,domblock312small(n),ierr)
        call MPI_Type_commit(domblock312small(n),ierr)

    enddo
    call MPI_Type_CONTIGUOUS(n2,MPI_REAL,columntype,ierr)
    call MPI_Type_commit(columntype,ierr)
    call MPI_Type_CONTIGUOUS(2*n2,MPI_REAL,columntype2,ierr)
    call MPI_Type_commit(columntype2,ierr)


end subroutine initializedomain
!**********************************************************************

subroutine dividedomain(n1,n2,ini)
    implicit none
    integer :: n1,n2
    integer :: ini(0:numtasks-1)
    real, dimension(n2big,n3big) :: varreadbig
    real, dimension(n1,n2) :: varread
    integer,  dimension(n2) :: ncells
    integer :: ntotal,ncount,i,j,n
    real :: nperpid

    write(6,*)'reading soil data',n1,n2,n2big,n3big
    !read in topo data to use as mask

    call sleep(10)
    if (PLATFORM_TO_RUN == 2) then
        open(11,file='/mnt/research/water/LHF_DATA/initialization/leafhydro_input_LE.dat',&
            form='unformatted',access='direct',recl=4*n1*n2)
    elseif (PLATFORM_TO_RUN == 1) then
        open(11,file='/egr/research-hydro/LHF_INI/initialization/leafhydro_input_LE.dat',&
            form='unformatted',access='direct',recl=4*n1*n2)
    elseif (PLATFORM_TO_RUN == 3) then
        open(11,file='/glade/p_old/work/shinsa11/DATA/CONUS_LHF_DATA/leafhydro_input_LE.dat',&
            form='unformatted',access='direct',recl=4*n1*n2)
    else
        stop "not supportted PLATFORM_TO_RUN"
    endif
    
    read(11,rec=2)((varreadbig(i,j),i=1,n2big),j=1,n3big)
    varread(1:n1,1:n2)=varreadbig(nw:ne,ns:nn)
    close(11)

    ntotal=count(varread>=0.1)
    nperpid=float(ntotal)/float(numtasks-2)
    write(6,*)'total number of land cells',ntotal


    ncells=count(varread>=0.1,1)

    ncount=0
    ini(1)=1
    n=2
    do j=1,n2
        ncount=ncount+ncells(j)
        if(ncount.ge.nint(float(n-1)*nperpid))then
            ini(n)=j-1
            n=n+1
        endif
        if(n.eq.numtasks-1)exit
    enddo

end subroutine dividedomain

!     ******************************************************************
subroutine SENDBORDERS(n2,js,je,wtd,reqsu,reqsd,reqru,reqrd)
    integer :: n2,js,je,reqsu,reqsd,reqru,reqrd
    real, dimension(n2,js:je):: wtd

    if(pid.eq.1)then
        call MPI_isend(wtd(1,je-1),1,columntype,2,200,MPI_COMM_WORLD,reqsu,ierr)
        call MPI_irecv(wtd(1,je),1,columntype,2,201,MPI_COMM_WORLD,reqru,ierr)

    elseif(pid.eq.numtasks-2)then
        call MPI_isend(wtd(1,js+1),1,columntype,pid-1,201,MPI_COMM_WORLD,reqsd,ierr)
        call MPI_irecv(wtd(1,js),1,columntype,pid-1,200,MPI_COMM_WORLD,reqrd,ierr)

    elseif(pid.gt.1.and.pid.lt.numtasks-2)then
        call MPI_isend(wtd(1,je-1),1,columntype,pid+1,200,MPI_COMM_WORLD,reqsu,ierr)
        call MPI_isend(wtd(1,js+1),1,columntype,pid-1,201,MPI_COMM_WORLD,reqsd,ierr)
        call MPI_irecv(wtd(1,js),1,columntype,pid-1,200,MPI_COMM_WORLD,reqrd,ierr)
        call MPI_irecv(wtd(1,je),1,columntype,pid+1,201,MPI_COMM_WORLD,reqru,ierr)

    endif

end subroutine sendborders

!     ******************************************************************
subroutine SENDBORDERSFLOOD(n2,js,je,wtd)
    integer :: n2,js,je
    real, dimension(2,n2,js:je):: wtd

    if(pid.eq.1)then
        call MPI_send(wtd(1,1,je-1),1,columntype2,2,400,MPI_COMM_WORLD,ierr)
        call MPI_recv(wtd(1,1,je),1,columntype2,2,401,MPI_COMM_WORLD,status,ierr)

    elseif(pid.eq.numtasks-2)then
        call MPI_send(wtd(1,1,js+1),1,columntype2,pid-1,401,MPI_COMM_WORLD,ierr)
        call MPI_recv(wtd(1,1,js),1,columntype2,pid-1,400,MPI_COMM_WORLD,status,ierr)

    elseif(pid.gt.1.and.pid.lt.numtasks-2)then
        call MPI_send(wtd(1,1,je-1),1,columntype2,pid+1,400,MPI_COMM_WORLD,ierr)
        call MPI_send(wtd(1,1,js+1),1,columntype2,pid-1,401,MPI_COMM_WORLD,ierr)
        call MPI_recv(wtd(1,1,js),1,columntype2,pid-1,400,MPI_COMM_WORLD,status,ierr)
        call MPI_recv(wtd(1,1,je),1,columntype2,pid+1,401,MPI_COMM_WORLD,status,ierr)

    endif

    end subroutine sendbordersflood

    !     ******************************************************************
    subroutine SENDFLOODBORDER(n2,js,je,var,varsouth,varnorth,reqsu,reqsd,reqru,reqrd)
    integer :: n2,js,je,reqsu,reqsd,reqru,reqrd
    real, dimension(n2,js:je):: var
    real, dimension(n2) :: varnorth,varsouth

    if(pid.eq.1)then
        call MPI_isend(var(1,je),1,columntype,2,300,MPI_COMM_WORLD,reqsu,ierr)
        call MPI_irecv(varnorth(1),1,columntype,2,301,MPI_COMM_WORLD,reqru,ierr)
        varsouth=0.

    elseif(pid.eq.numtasks-2)then
        call MPI_isend(var(1,js),1,columntype,pid-1,301,MPI_COMM_WORLD,reqsd,ierr)
        call MPI_irecv(varsouth(1),1,columntype,pid-1,300,MPI_COMM_WORLD,reqrd,ierr)
        varnorth=0.

    elseif(pid.gt.1.and.pid.lt.numtasks-2)then
        call MPI_isend(var(1,je),1,columntype,pid+1,300,MPI_COMM_WORLD,reqsu,ierr)
        call MPI_isend(var(1,js),1,columntype,pid-1,301,MPI_COMM_WORLD,reqsd,ierr)
        call MPI_irecv(varsouth(1),1,columntype,pid-1,300,MPI_COMM_WORLD,reqrd,ierr)
        call MPI_irecv(varnorth(1),1,columntype,pid+1,301,MPI_COMM_WORLD,reqru,ierr)

    endif

end subroutine sendfloodborder

END MODULE MODULE_PARALLEL
