program main

!*****************************************************************************80
!! Calculate Lyapunove exponent of Lorenz 63 model
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 July 2022
!
!  Author:
!
!    Mao Ouyang
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  integer, parameter :: n = 3
  integer, parameter :: m = 8000000
  integer :: i, j

  external l63_dydt
  real ( kind = rk ) :: y0(n)
  real ( kind = rk ) :: dt = 0.01
  real ( kind = rk ) :: tspan(2)
  real ( kind = rk ) :: t(m + 1)
  real ( kind = rk ) :: y(m + 1, n)

  character ( len = 255 ) data_filename
  integer data_unit

  call timestamp ( )
  write ( *, '(a)' ) 'Nature run of Lorenz 63 in CSE (MS2022)'
  y0(1) = 8.20747939
  y0(2) = 10.0860429
  y0(3) = 23.0860429
  write ( *, '(a)' ) 'Initial state variables'
  write(*, '(3(F12.8))') y0

  tspan(1) = 0.0
  tspan(2) = 80000
  call rk4 (l63_dydt, tspan, y0, m, n, t, y)

  write (*, *) 'Save the data'

  call get_unit ( data_unit )

  data_filename = 'Lorenz63_NR.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do i = 1, m + 1
    write ( data_unit, '(4(f15.8))' ) t(i), y(i, :)
  end do

  close ( unit = data_unit )

  write ( *, '(a)' ) 'Normal end of execution.'
  call timestamp ( )

end

subroutine l63_dydt( t, y, value )

!*****************************************************************************80
!
!! l63_dydt returns the right hand side of the Lorenz 63 model.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2022
!
!  Author:
!
!    Mao Ouyang
!
!  Input:
!
!    real ( kind = rk ) t, y(:): the arguments.
!
!  Output:
!
!    real ( kind = rk ) value(:): the value of the derivative.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 ), n = 3

  real ( kind = rk ) t
  real ( kind = rk ) value(n)
  real ( kind = rk ) y(n)
  real ( kind = rk ) alpha, rho, beta

  call r8_fake_use ( t )

  alpha = 10.0
  rho = 28.0
  beta = 8.0 / 3.0

  value(1) = alpha * (y(2) - y(1))
  value(2) = -y(1) * y(3) + rho * y(1) - y(2)
  value(3) = y(1) * y(2) - beta * y(3)

  return
end
