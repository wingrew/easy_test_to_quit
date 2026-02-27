!-----------------------------------------------------------------------------------------------------------------
! three dimensional interpolation for vertex center grid structure  
  subroutine global_interp(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA

!~~~~~~> Other parameters:

  integer :: j,m,imin,jmin,kmin
  integer,dimension(3) :: cxB,cxT,cxI,cmin,cmax
  real*8,dimension(3) :: cx
  real*8, dimension(1:ORDN) :: x1a
  real*8, dimension(1:ORDN,1:ORDN,1:ORDN) :: ya
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: dX,dY,dZ,ddy
  real*8, parameter :: ONE=1.d0
  logical::decide3d

  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  dX = X(imin+1)-X(imin)
  dY = Y(jmin+1)-Y(jmin)
  dZ = Z(kmin+1)-Z(kmin)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  cxI(1) = idint((x1-X(1))/dX+0.4)+1
  cxI(2) = idint((y1-Y(1))/dY+0.4)+1
  cxI(3) = idint((z1-Z(1))/dZ+0.4)+1

  cxB = cxI - ORDN/2+1
  cxT = cxB + ORDN - 1
       
  cmin = 1
  cmax = ex
  if(Symmetry == OCTANT  .and.dabs(X(1))<dX) cmin(1) = -ORDN/2+2
  if(Symmetry == OCTANT  .and.dabs(Y(1))<dY) cmin(2) = -ORDN/2+2
  if(Symmetry /= NO_SYMM .and.dabs(Z(1))<dZ) cmin(3) = -ORDN/2+2
  do m =1,3
   if(cxB(m) < cmin(m))then
      cxB(m) = cmin(m)
      cxT(m) = cxB(m) + ORDN - 1
   endif
   if(cxT(m) > cmax(m))then
      cxT(m) = cmax(m)
      cxB(m) = cxT(m) + 1 - ORDN
   endif
 enddo
 if(cxB(1)>0)then
  cx(1) = (x1 - X(cxB(1)))/dX
 else
  cx(1) = (x1 + X(2-cxB(1)))/dX
 endif
 if(cxB(2)>0)then
  cx(2) = (y1 - Y(cxB(2)))/dY
 else
  cx(2) = (y1 + Y(2-cxB(2)))/dY
 endif
 if(cxB(3)>0)then
  cx(3) = (z1 - Z(cxB(3)))/dZ
 else
  cx(3) = (z1 + Z(2-cxB(3)))/dZ
 endif

  if(decide3d(ex,f,f,cxB,cxT,SoA,ya,ORDN,Symmetry))then
     write(*,*)"global_interp position: ",x1,y1,z1
     write(*,*)"data range: ",X(1),X(ex(1)),Y(1),Y(ex(2)),Z(1),Z(ex(3))
     stop
  endif
  call polin3(x1a,x1a,x1a,ya,cx(1),cx(2),cx(3),f_int,ddy,ORDN)

  return

  end subroutine global_interp