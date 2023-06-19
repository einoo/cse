program main

!*****************************************************************************80
!
!! Generate random normal distribution with mean and standard deviations
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

  integer data_unit, i

  integer, parameter :: m = 8000001
  integer, parameter :: n = 3

  real ( kind = rk ) mu
  real ( kind = rk ) r(m,n)
  real ( kind = rk ) y(m,n+1)
  real ( kind = rk ) sigma

  write ( *, '(a)' ) 'R8MAT_NORMAL_AB returns a matrix of Normal AB values.'

  mu = 0.0D+00
  sigma = SQRT(2.0D+00)

  write ( *, '(a,g14.6)' ) '  MU = ', mu
  write ( *, '(a,g14.6)' ) '  SIGMA = ', sigma

  call r8mat_normal_ab ( m, n, mu, sigma, r )

  write ( *, * ) 'Read the true data'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = '../Data/Lorenz63_NR.txt', status = 'old', &
         access = 'sequential', form='formatted', action='read')
  do i = 1, m
    read( data_unit, '(4(f15.8))') y(i, :)
  end do
  close (unit = data_unit)

  write ( *, * ) 'Save the observation data'
  call l96_observation_save ( m, n, y, r)

end

subroutine l96_observation_save ( m, n, y, r )

!*****************************************************************************80
!
!! Save the lorenz 96 observations (true + normal distribution).
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
!  Input:
!
!    integer m: the number of steps to take.
!
!    integer n: the number of variables.
!
!    real ( kind = rk ) y(m, n+1), r(m, n): the true and normal distributions.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  character ( len = 255 ) data_filename
  integer data_unit
  integer i, j
  real ( kind = rk ) y(m, n+1)
  real ( kind = rk ) r(m, n)
!
!  Create the data file.
!
  call get_unit ( data_unit )

  write (*, *) 'Save the noise for normal distribution verification'
  open ( unit = data_unit, file = 'noise.txt', status = 'replace' )
  do i = 1, m
    write (*, *) i
    write ( data_unit, '(4(f15.8))' ) r(i, :)
  end do
  close ( unit = data_unit )

  data_filename = 'Lorenz63_OBS.txt'

  open ( unit = data_unit, file = data_filename, status = 'replace' )

  do j = 2, n + 1
    do i = 1, m
      r (i, j-1) = y(i, j) + r(i, j-1)
      write(*, *) i, j
    end do
  end do


  do i = 1, m
    write (*, *) i
    write ( data_unit, '(4(f15.8))' ) y(i, 1), r(i, :)
  end do

  close ( unit = data_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  data stored in "' &
    // trim ( data_filename ) // '".'
  return
end
