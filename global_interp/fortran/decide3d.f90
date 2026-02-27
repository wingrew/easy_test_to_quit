  function decide3d(ex,f,fpi,cxB,cxT,SoA,ya,ORDN,Symmetry)  result(gont)
  implicit none

  integer,                                 intent(in) :: ORDN,Symmetry
  integer,dimension(1:3)                 , intent(in) :: ex,cxB,cxT
  real*8, dimension(1:3)                 , intent(in) :: SoA
  real*8, dimension(ex(1),ex(2),ex(3))   , intent(in) :: f,fpi
  real*8, dimension(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3)), intent(out):: ya
  logical::gont

  integer,dimension(1:3) :: fmin1,fmin2,fmax1,fmax2
  integer::i,j,k,m

  gont=.false.
  do m=1,3
! check cxB and cxT are NaN or not  
    if(.not.(iabs(cxB(m)).ge.0)) gont=.true.
    if(.not.(iabs(cxT(m)).ge.0)) gont=.true.
    fmin1(m) = max(1,cxB(m))
    fmax1(m) = cxT(m)
    fmin2(m) = cxB(m)
    fmax2(m) = min(0,cxT(m))
    if((fmin1(m).le.fmax1(m)).and.(  fmin1(m)<1.or.  fmax1(m)>ex(m)))gont=.true.
    if((fmin2(m).le.fmax2(m)).and.(2-fmax2(m)<1.or.2-fmin2(m)>ex(m)))gont=.true.
  enddo
!sanity check
  if(gont)then
          write(*,*)"error in decide3d"
          write(*,*)((fmin1.le.fmax1).and.(  fmin1<1.or.  fmax1>ex))
          write(*,*)((fmin2.le.fmax2).and.(2-fmax2<1.or.2-fmin2>ex))
          write(*,*)"cxB, cxT and data shape:"
          write(*,*)cxB,cxT,ex
          write(*,*)"resulted fmin1, fmax1 and fmin2, fmax2:"
          write(*,*)fmin1,fmax1,fmin2,fmax2
  else

  do k=fmin1(3),fmax1(3)
     do j=fmin1(2),fmax1(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,k)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,j,k)*SoA(1)
        enddo
     enddo
     do j=fmin2(2),fmax2(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,2-j,k)*SoA(2)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,2-j,k)*SoA(1)*SoA(2)
        enddo
     enddo
  enddo
           
  do k=fmin2(3),fmax2(3)
     do j=fmin1(2),fmax1(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,2-k)*SoA(3)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,j,2-k)*SoA(1)*SoA(3)
        enddo
     enddo
     do j=fmin2(2),fmax2(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,2-j,2-k)*SoA(2)*SoA(3)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,2-j,2-k)*SoA(1)*SoA(2)*SoA(3)
        enddo
     enddo
  enddo

  endif

  end function decide3d