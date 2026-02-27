program benchmark_fortran
  implicit none
  integer :: nquery
  character(len=512) :: arg
  character(len=512) :: out_path
  integer, parameter :: ex1=24, ex2=20, ex3=18, ordn=4
  integer :: ex(3), symmetry
  real*8 :: soa(3)
  real*8 :: X(ex1), Y(ex2), Z(ex3)
  real*8 :: f(ex1,ex2,ex3)
  real*8 :: x_max, y_max, z_max
  integer :: i,j,k,q
  integer :: unit_id, ios
  integer :: count0, count1, rate
  real*8 :: elapsed, out, checksum
  integer(kind=8) :: state
  real*8 :: rx, ry, rz, xq, yq, zq

  call get_command_argument(1, arg)
  read(arg,*) nquery
  call get_command_argument(2, out_path)

  ex = (/ex1, ex2, ex3/)
  symmetry = 2
  soa = (/1.d0, -1.d0, 1.d0/)

  do i=1,ex1
    X(i) = dble(i-1) * 0.1d0
  end do
  do j=1,ex2
    Y(j) = dble(j-1) * 0.2d0
  end do
  do k=1,ex3
    Z(k) = dble(k-1) * 0.15d0
  end do

  do k=1,ex3
    do j=1,ex2
      do i=1,ex1
        f(i,j,k) = dsin(X(i)) + 0.5d0 * dcos(2.d0 * Y(j)) + Z(k)*Z(k) + 0.1d0*X(i)*Y(j)*Z(k)
      end do
    end do
  end do

  x_max = X(ex1)
  y_max = Y(ex2)
  z_max = Z(ex3)

  open(newunit=unit_id, file=trim(out_path), status='replace', action='write', iostat=ios)
  if (ios /= 0) then
    write(*,*) 'failed to open output file: ', trim(out_path)
    stop 2
  end if

  checksum = 0.d0
  state = 123456789_8

  call system_clock(count0, rate)
  do q=1,nquery
    state = mod(1664525_8 * state + 1013904223_8, 4294967296_8)
    rx = dble(state) / 4294967295.d0
    state = mod(1664525_8 * state + 1013904223_8, 4294967296_8)
    ry = dble(state) / 4294967295.d0
    state = mod(1664525_8 * state + 1013904223_8, 4294967296_8)
    rz = dble(state) / 4294967295.d0

    xq = rx * x_max
    yq = ry * y_max
    zq = rz * z_max

    call global_interp(ex, X, Y, Z, f, out, xq, yq, zq, ordn, soa, symmetry)
    checksum = checksum + out
    write(unit_id, '(ES26.17)') out
  end do
  call system_clock(count1)

  close(unit_id)

  elapsed = dble(count1 - count0) / dble(rate)
  write(*,'(A,F0.9)') 'elapsed_sec=', elapsed
  write(*,'(A,ES26.17)') 'checksum=', checksum
end program benchmark_fortran
