program main

!*****************************************************************************80
!
!! Ensemble Kalman Filter (Perturb Observation)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2022
!
!  Author:
!
!    Mao Ouyang
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  integer, parameter :: n = 3  ! state variables
  integer, parameter :: m = 3  ! ensemble size
  integer, parameter :: ta = 8  ! time step for observation assimilation
  integer, parameter :: nt = 20800  ! data assimilation steps
  real ( kind = rk ) :: dt = 0.01 ! time step

  external l63_dydt
  real ( kind = rk ) :: tspan(2)
  real ( kind = rk ) :: t(ta + 1)
  real ( kind = rk ) :: y(ta + 1, n)
  real ( kind = rk ) :: rmsea(nt)
  real ( kind = rk ) :: rmseo(nt)
  integer :: i, j

  integer data_unit
  character ( len = 255 ) data_filename

  real ( kind = rk ) :: x0(m, n+1)  ! initial ensemble
  real ( kind = rk ) :: delta ! inflation variance
  real ( kind = rk ) :: h(n, n) ! observation matrix
  real ( kind = rk ) :: r(n, n) ! observation error covariance matrix
  real ( kind = rk ) :: xtru(nt, n+1) ! true data
  real ( kind = rk ) :: xobs(nt, n+1) ! observation data
  real ( kind = rk ) :: mu, sigma
  real ( kind = rk ) :: xf(m, n)  ! background/forecast state variables
  real ( kind = rk ) :: xa(nt, n)  ! analysis/predict state variables
  real ( kind = rk ) :: xas(nt, m, n)  ! analysis/predict state variables (save data)
  real ( kind = rk ) :: xfs(nt, m, n)  ! forecast ensembles (save data)
  real ( kind = rk ) :: xda(m, n)  ! analysis/predict ensembles

  call timestamp ( )

  write (*, *) 'Load the true and observation data'
  xtru(:, :) = 0.0
  call get_unit ( data_unit )
  open ( unit = data_unit, file = '../Data/Lorenz63_NR.txt', status = 'old', &
         access = 'sequential', form='formatted', action='read')
  do i = 1, nt
    read( data_unit, '(4(f15.8))') xtru(i, :)
  end do
  close (unit = data_unit)

  call get_unit ( data_unit )
  open ( unit = data_unit, file = '../Data/Lorenz63_OBS.txt', status = 'old', &
         access = 'sequential', form='formatted', action='read')
  do i = 1, nt
    read( data_unit, '(4(f15.8))') xobs(i, :)
  end do
  close (unit = data_unit)

  write (*, *) 'Observation matrix and error'
  h = 0
  forall(i=1:n) h(i, i) = 1
  r = 0
  forall(i=1:n) r(i, i) = 2

  write (*, *) 'Read the initial ensembles'
  do i = 1, m
    call get_unit ( data_unit )
    write(data_filename, '(a, i0.6, a)') '../Data/Init/', i, '/x0'
    open ( unit = data_unit, file = data_filename, status = 'old', &
           access = 'sequential', form='formatted', action='read')
    read( data_unit, '(4(f15.8))') x0(i, :)
    close (unit = data_unit)
  end do

  write (*, *) 'Ensemble Transform Kalman filter (ETKF) MS2022 '
  xf(:, :) = x0(:, 2:n+1)
  delta = 0.04
  xda(:, :) = 0.0
  do i = 1, nt, ta
    write (*, *) i, m, n
    xfs(i, :, :) = xf(:, :)
    write (*, *) 'True data'
    write(*, '(3(f15.8))') xtru(i, 2:n+1)
    write (*, *) 'Observation data'
    write(*, '(3(f15.8))') xobs(i, 2:n+1)
    call etkf (xf, m, n, delta, h, r, xobs(i, 2:), xda)
    xas(i, :, :) = xda(:, :)
    ! calculate the analysis ensemble mean
    do j = 1, n
      xa(i, j) = SUM(xda(:, j)) / m
    end do
    write (*, *) 'Analysis/update mean'
    write (*, '(3(f15.8))') xa(i, :)
    write (*, *) 'RMSE'
    write (*, *) SQRT(SUM((xa(i, :) - xtru(i, 2:4)) ** 2) / n)
    rmsea(i) = SQRT(SUM((xa(i, :) - xtru(i, 2:4)) ** 2) / n)
    rmseo(i) = SQRT(SUM((xobs(i, 2:4) - xtru(i, 2:4)) ** 2) / n)
    ! free forecast for ta steps
    do j = 1, m
      tspan(1) = i * dt
      tspan(2) = (i + 8) * dt
      call rk4 (l63_dydt, tspan, xda(j, :), ta, n, t, y)
      xf(j, :) = y(ta+1, :)
    end do
  end do

  ! save the data
  call get_unit ( data_unit )
  open ( unit = data_unit, file = '../Data/rmse.txt', status='replace')
  do i = 1, nt, ta
    write (data_unit, '(2(f15.8))') rmseo(i), rmsea(i)
  end do
  close (unit = data_unit)

  call get_unit ( data_unit )
  open ( unit = data_unit, file = '../Data/xa_full.txt', status='replace')
  do i = 1, nt, ta
    write (data_unit, '(3000(f15.8))') xas(i, :, :)
  end do
  close (unit = data_unit)

  call get_unit ( data_unit )
  open ( unit = data_unit, file = '../Data/xf_full.txt', status='replace')
  do i = 1, nt, ta
    write (data_unit, '(3000(f15.8))') xfs(i, :, :)
  end do
  close (unit = data_unit)

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
