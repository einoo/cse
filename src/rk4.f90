subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer ios
  integer iunit
  logical lopen
 
  iunit = 0
 
  do i = 1, 99
 
    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )
 
      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if
 
  end do

  return
end
subroutine r8_fake_use ( x )

!*****************************************************************************80
!
!! r8_fake_use pretends to use a variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  r8_fake_use: variable is NAN.'
  end if

  return
end
subroutine rk4 ( dydt, tspan, y0, n, m, t, y )

!*****************************************************************************80
!
!! rk4 approximates an ODE using a Runge-Kutta fourth order method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    function handle dydt: a function that evaluates the right hand side.
!
!    real ( kind = rk ) tspan(2): contains the initial and final times.
!
!    real ( kind = rk ) y0(m): the initial condition.
!
!    integer n: the number of steps to take.
!
!    integer m: the number of variables.
!
!  Output:
!
!    real ( kind = rk ) t(n+1), y(n+1,m): the times and solution values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) dt
  external dydt
  real ( kind = rk ) f1(m)
  real ( kind = rk ) f2(m)
  real ( kind = rk ) f3(m)
  real ( kind = rk ) f4(m)
  integer i
  real ( kind = rk ) t(n+1)
  real ( kind = rk ) tspan(2)
  real ( kind = rk ) y(n+1,m)
  real ( kind = rk ) y0(m)

  dt = ( tspan(2) - tspan(1) ) / real ( n, kind = rk )

  t(1) = tspan(1)
  y(1,1:m) = y0(1:m)

  do i = 1, n

    call dydt ( t(i),            y(i,1:m),                 f1 )
    call dydt ( t(i) + dt / 2.0, y(i,1:m) + dt * f1 / 2.0, f2 )
    call dydt ( t(i) + dt / 2.0, y(i,1:m) + dt * f2 / 2.0, f3 )
    call dydt ( t(i) + dt,       y(i,1:m) + dt * f3,       f4 )

    t(i+1) = t(i) + dt
    y(i+1,1:m) = y(i,1:m) + dt * ( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine rk4_tml ( dydt, tspan, y0, n, m, a, t, y )

!*****************************************************************************80
!
!! rk4 approximates an ODE using a Runge-Kutta fourth order method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2020
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    function handle dydt: a function that evaluates the right hand side.
!
!    real ( kind = rk ) tspan(2): contains the initial and final times.
!
!    real ( kind = rk ) y0(m): the initial condition.
!
!    real ( kind = rk ) a(m): the state variables for calculating the PDF.
!
!    integer n: the number of steps to take.
!
!    integer m: the number of variables.
!
!  Output:
!
!    real ( kind = rk ) t(n+1), y(n+1,m): the times and solution values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) dt
  external dydt
  real ( kind = rk ) f1(m, m)
  real ( kind = rk ) f2(m, m)
  real ( kind = rk ) f3(m, m)
  real ( kind = rk ) f4(m, m)
  real ( kind = rk ) a(m)
  integer i, j
  real ( kind = rk ) t(n+1)
  real ( kind = rk ) tspan(2)
  real ( kind = rk ) y(n+1,m, m)
  real ( kind = rk ) y0(m, m)

  dt = ( tspan(2) - tspan(1) ) / real ( n, kind = rk )

  t(1) = tspan(1)
  y(1,:, :) = y0

  do i = 1, n

    call dydt ( t(i),            y(i, :, :),                 a, f1 )
    call dydt ( t(i) + dt / 2.0, y(i, :, :) + dt * f1 / 2.0, a, f2 )
    call dydt ( t(i) + dt / 2.0, y(i, :, :) + dt * f2 / 2.0, a, f3 )
    call dydt ( t(i) + dt,       y(i, :, :) + dt * f3,       a, f4 )

    t(i+1) = t(i) + dt
    y(i+1, :, :) = y(i, :, :) + dt * ( f1(:, :) + 2.0 * f2(:, :) + 2.0 * f3(:, :) + f4(:, :) ) / 6.0

  end do

  return
end
