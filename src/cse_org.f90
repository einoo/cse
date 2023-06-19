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

  integer, parameter :: rk = kind(1.0D+00)
  integer, parameter :: n = 3  ! state variables
  integer, parameter :: m = 3  ! ensemble size
  integer, parameter :: ta = 8  ! time step for observation assimilation
  integer, parameter :: nt = 20800  ! data assimilation steps
  ! integer, parameter :: nt = 10000 ! data assimilation steps
  real(kind=rk), parameter :: t0 = 75.1  ! time period
  integer, parameter :: ms = NINT(t0 * 4) ! time period for free forecast
  real(kind=rk), parameter :: dt = 0.01 ! time step

  external l63_dydt
  real(kind=rk) :: tspan(2)
  real(kind=rk) :: y0(n)
  real(kind=rk) :: t(ms + 1), tt(2), tc(nt + 1)
  real(kind=rk) :: y(ms + 1, n), yt(2, n), yc(nt + 1, n)
  real(kind=rk) :: ys(m, ms + 1, n) ! saving the free forecast for control
  real(kind=rk) :: yse(m, ms + ta + 1, n) ! saving the free forecast for extended control
  real(kind=rk) :: rmsea(nt)
  real(kind=rk) :: rmseo(nt)
  ! integer :: i, j, k, info, ts=(406-1) * 8 + 1 ! initial time step for CSE
  integer :: i, j, k, info, ts=(1-1) * 8 + 1 ! initial time step for CSE

  real(kind=rk), allocatable :: te(:), ye(:, :)

  integer data_unit, pens, nrs, srs, backn, srs2
  character(len=255) data_filename

  real(kind=rk) :: x0(n*m)
  real(kind=rk) :: delta ! inflation variance
  real(kind=rk) :: h(n, n) ! observation matrix
  real(kind=rk) :: r(n, n) ! observation error covariance matrix
  real(kind=rk) :: xtru(n + 1) ! true data at the initial state
  real(kind=rk) :: xobs(n + 1) ! observation data at the initial state
  real(kind=rk) :: mu, sigma
  real(kind=rk) :: xf(m, n)  ! background/forecast state variables
  real(kind=rk) :: xa(nt, n)  ! analysis/predict state variables
  real(kind=rk) :: xas(nt, m, n)  ! analysis/predict state variables (save data)
  real(kind=rk) :: xfs(nt, m, n)  ! forecast ensembles (save data)
  real(kind=rk) :: xda(m, n)  ! analysis/predict ensembles
  real(kind=rk) :: obs_noise(nt, n)  ! observation noise

  real(kind=rk) :: d = 0.05 ! Euclidean norm
  real(kind=rk) :: dc(ta - 1, n) ! perturbation added to the state variables


  call timestamp()

  write (*, *) 'Reproduce the Control System Experiment of Myoshi and Sun 2022'
  write (*, *) ' '
  write (*, *) 'Load the true and observation data at the initial step'
  call get_unit(data_unit)
  open (unit=data_unit, file='xtrue', status='old', &
        access='sequential', form='formatted', action='read')
  read (data_unit, '(4(f15.8))') xtru
  close (unit=data_unit)
  write (*, '(4(f15.8))') xtru

  call get_unit(data_unit)
  open (unit=data_unit, file='xobs', status='old', &
        access='sequential', form='formatted', action='read')
  read (data_unit, '(4(f15.8))') xobs
  close (unit=data_unit)
  write (*, '(4(f15.8))') xobs

  write (*, *) 'Load the observation noise'
  call get_unit(data_unit)
  open (unit=data_unit, file='../Data/noise.txt', status='old', &
        access='sequential', form='formatted', action='read')
  do i = 1, nt
    read (data_unit, '(3(f15.8))') obs_noise(i, :)
  end do
  close (unit=data_unit)
  write (*, '(4(f15.8))') xobs(1), obs_noise(ts, :)

  write (*, *) 'Observation location and error covariance matrix'
  h = 0
  forall (i=1:n) h(i, i) = 1
  r = 0
  forall (i=1:n) r(i, i) = 2

  write (*, *) 'Read the initial ensembles'
  call get_unit(data_unit)
  open (unit=data_unit, file='xf', status='old', &
        access='sequential', form='formatted', action='read')
  read (data_unit, '(9(f15.8))') x0
  close (unit=data_unit)
  xf = RESHAPE(x0, (/m, n/))
  write (*, '(3(f15.8))') xf

  write (*, *) 'Ensemble Transform Kalman filter (ETKF) MS2022 '
  delta = 0.04
  xda(:, :) = 0.0
  do i = ts, nt, ta
    write (*, *) i, m, n
    xfs(i, :, :) = xf(:, :)
    write (*, *) 'True data'
    write (*, '(3(f15.8))') xtru(2:n + 1)
    write (*, *) 'Observation data'
    write (*, '(3(f15.8))') xobs(2:n + 1)
    call etkf(xf, m, n, delta, h, r, xobs(2:n + 1), xda)
    xas(i, :, :) = xda(:, :)
    ! calculate the analysis ensemble mean
    do j = 1, n
      xa(i, j) = SUM(xda(:, j))/m
    end do
    write (*, *) 'Analysis/update mean'
    write (*, '(3(f15.8))') xa(i, :)
    write (*, *) 'RMSE'
    write (*, *) SQRT(SUM((xa(i, :) - xtru(2:n + 1))**2)/n)
    rmsea(i) = SQRT(SUM((xa(i, :) - xtru(2:n + 1))**2)/n)
    rmseo(i) = SQRT(SUM((xobs(2:n + 1) - xtru(2:n + 1))**2)/n)
    ! free forecast for T steps
    write (*, *) 'Step 2 Free forecast for T steps'
    write (*, '(a)') 'Initial state variables'
    do j = 1, m
      tspan(1) = 0.0
      tspan(2) = 3.0
      write (*, '(3(f15.8))') xda(j, :)
      call rk4(l63_dydt, tspan, xda(j, :), ms, n, t, y)
      ys(j, :, :) = y(:, :)
    end do
    ! determine whether activating control or not
    call cse_check (ys(:, :, :), m, n, ms, d, ta, info, dc, srs, nrs)
    write (*, '(i3)') info
    write (*, '(3(f15.8))') dc
    if (allocated(te)) deallocate(te)
    if (allocated(ye)) deallocate(ye)
    pens = i - ta
    do while (info .eq. -1)
      write (*, *) 'COME TO HERE'
      backn = 1 ! the number of loop in finding former ensembles
      write (*, '((i10), (a), (i10))') i, '  ALL ENSEMBLES SHOW REGIME SHIFT', pens
      write (*, '((a), (i1), (i1), (i1))') 'Order of ensemble shows regime shift is: ', srs, nrs, 6-srs-nrs
      srs2 = nrs ! find a second ensemble member showed regime shift
      ! Use the previous state variables to find a ensemble without regime shift
      allocate (te(ms + ta*backn + 1))
      allocate (ye(ms + ta*backn + 1, n))
      do j = 1, m
        tspan(1) = 0.0
        tspan(2) = 3.0 + dt*backn*ta
        call rk4(l63_dydt, tspan, xas(pens, j, :), ms + ta*backn, n, te, ye)
        if (all(ye((ta*backn):(ms + ta*backn + 1), 1) .ge. 0)) then
          info = 1
          nrs = j
          exit
        else
          info = -1
          nrs = 0
        end if
      end do
      write (*, *) pens, info, nrs, srs
      if  (srs .eq. nrs) srs = srs2
      ! write (*, *) srs, nrs, ys(srs, ms+backn*ta+1, 1), ys(nrs, ms+backn*ta+1, 1)
      do j = 1, ta - 1
        dc(j, :) = ys(nrs, j, :) - ys(srs, j, :)
        dc(j, :) = dc(j, :)/NORM2(dc(j, :))*d
      end do

      if (info .eq. -1) then
        pens = pens - ta
        backn = backn + 1
      end if
      deallocate(te, ye)
      if (info .ne. -1) then
        write (*, '(a)') 'Successfully find the ensemble not showing regime shift'
        write (*, *) i, info, nrs, srs
      end if
    end do

    ! conduct the control
    y0 = xtru(2:4)
    yc(i, :) = xtru(2:4)
    tc(i) = i * dt
    do j = 1, ta-1
      tspan(1) = 0.0
      tspan(2) = dt
      call rk4(l63_dydt, tspan, y0, 1, n, tt, yt)
      yc(i + j, :) = yt(2, :)
      tc(i + j) = (i + j)*dt
      y0 = yt(2, :) + dc(j, :)
      if (yc(i + j, 1) .le. 0) then
        write (*, *) yc(i + j, 1), i, j
        ! stop
      end if
    end do
    ! run another free forecast to obtain the true data for the next DA cycle
    ! observation is generated based on the new true data with noise
    tspan(1) = 0.0
    tspan(2) = dt
    call rk4(l63_dydt, tspan, y0, 1, n, tt, yt)
    write (*, *) 'Controlled nature run'
    do j = 0, ta -1
      write (*, '(3(f15.8))') yc(i + j, :)
    end do

    ! assign the values to the next step
    xf = ys(:, 9, :)

    xtru(2:4) = yt(2, :)
    xobs(2:4) = xtru(2:4) + obs_noise(i + ta, :)
    write (*, *) 'True'
    write (*, '(3(f15.8))') xtru(2:4)
    write (*, *) 'OBS'
    write (*, '(3(f15.8))') xobs(2:4)
    write (*, *) 'xf'
    write (*, '(3(f15.8))') xf
  end do

  ! save the data
  call get_unit(data_unit)
  open (unit=data_unit, file='cse_org.txt', status='replace')
  do i = ts, nt
    write (data_unit, '(4(f15.8))') tc(i), yc(i, :)
  end do
  close (unit=data_unit)

end

subroutine l63_dydt(t, y, value)

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

  integer, parameter :: rk = kind(1.0D+00), n = 3

  real(kind=rk) t
  real(kind=rk) value(n)
  real(kind=rk) y(n)
  real(kind=rk) alpha, rho, beta

  call r8_fake_use(t)

  alpha = 10.0
  rho = 28.0
  beta = 8.0/3.0

  value(1) = alpha*(y(2) - y(1))
  value(2) = -y(1)*y(3) + rho*y(1) - y(2)
  value(3) = y(1)*y(2) - beta*y(3)

  return
end
