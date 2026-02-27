!------------------------------------------------------------------------------
!
! interpolation in 3 dimensions, follow zyx order
!
!------------------------------------------------------------------------------
  subroutine polin3(x1a,x2a,x3a,ya,x1,x2,x3,y,dy,ordn)
  implicit none

  integer,intent(in) :: ordn
  real*8, dimension(1:ordn), intent(in) :: x1a,x2a,x3a
  real*8, dimension(1:ordn,1:ordn,1:ordn), intent(in) :: ya
  real*8, intent(in) :: x1,x2,x3
  real*8, intent(out) :: y,dy

  integer  :: j, k
  real*8, dimension(ordn,ordn) :: yatmp
  real*8, dimension(ordn) :: ymtmp
  real*8 :: dy_temp

  do k=1,ordn
    do j=1,ordn
      call polint(x1a, ya(:,j,k), x1, yatmp(j,k), dy_temp, ordn)
    end do
  end do
  do k=1,ordn
    call polint(x2a, yatmp(:,k), x2, ymtmp(k), dy_temp, ordn)
  end do
  call polint(x3a, ymtmp, x3, y, dy, ordn)


  return
  end subroutine polin3