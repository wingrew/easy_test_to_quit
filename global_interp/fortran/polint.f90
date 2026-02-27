  subroutine polint(xa, ya, x, y, dy, ordn)
  implicit none

  integer, intent(in) :: ordn
  real*8, dimension(ordn), intent(in) :: xa, ya
  real*8, intent(in) :: x
  real*8, intent(out) :: y, dy

  integer :: i, m, ns, n_m
  real*8, dimension(ordn) :: c, d, ho
  real*8 :: dif, dift, hp, h, den_val

  c = ya
  d = ya
  ho = xa - x

  ns = 1
  dif = abs(x - xa(1))

  do i = 2, ordn
    dift = abs(x - xa(i))
    if (dift < dif) then
      ns = i
      dif = dift
    end if
  end do

  y = ya(ns)
  ns = ns - 1

  do m = 1, ordn - 1
    n_m = ordn - m
    do i = 1, n_m
      hp = ho(i)
      h  = ho(i+m)
      den_val = hp - h

      if (den_val == 0.0d0) then
        write(*,*) 'failure in polint for point',x
        write(*,*) 'with input points: ',xa
        stop
      end if

      den_val = (c(i+1) - d(i)) / den_val

      d(i) = h * den_val
      c(i) = hp * den_val
    end do

    if (2 * ns < n_m) then
      dy = c(ns + 1)
    else
      dy = d(ns)
      ns = ns - 1
    end if
    y = y + dy
  end do

  return
  end subroutine polint

