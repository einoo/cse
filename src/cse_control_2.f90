subroutine cse_check(ys, m, n, ms, d, ta, info, dc, srs, nrs)

!*****************************************************************************80
!
!! check free forecast to check whether to perform control
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2022
!
!  Author:
!
!    Mao Ouyang
!
!  Input:
!
!    real ( kind = rk ) ys(ms+1, m, n): forecast matrix
!
!    integer n, m, ms, d, ta: the number of variables, ensembles, forecast time steps, norm, da time step
!
!  Output:
!
!    real ( kind = rk ) dc(n): control perturbation, the orthogonal matrix and L2 norm.
!
!    integer info : output information
!
  implicit none

  integer, parameter :: rk = kind(1.0D+00)

  integer :: n, m, ms, ta, info, i, j, srs, nrs

  real(kind=rk) :: ys(m, ms + 1, n)
  real(kind=rk) :: dc(ta - 1, n)
  real(kind=rk) :: gd(n)
  real(kind=rk) :: d

  if (all(ys(:, :, 1) .gt. 0)) then
      info = 0
      dc = 0
      srs = 0
      nrs = 0
  else
  do i = 1, ms
    gd = ys(:, i, 1) - ys(:, i-1, 1)
    if ((minval(ys(:, i, 1)) .lt. 0) &
        .and. (maxval(ys(:, i, 1)) .gt. 0) &
        .and. (minval(gd) .le. 0) &
        .and. (maxval(gd) .gt. -0.01)) then
      info = 1
      srs = minloc(ys(:, i, 1), DIM=1)
      nrs = maxloc(ys(:, i, 1), DIM=1)
      write (*, *) srs, nrs, ys(srs, i, 1), ys(nrs, i, 1)
      do j = 1, ta - 1
        dc(j, :) = ys(nrs, j, :) - ys(srs, j, :)
        dc(j, :) = dc(j, :)/NORM2(dc(j, :))*d
      end do
      exit
    else
      info = -1
      dc = 0.0
      srs = minloc(ys(:, ms+1, 1), DIM=1)
      nrs = 6 - srs - maxloc(ys(:, ms+1, 1), DIM=1)
    end if
  end do
  end if

  return
end
