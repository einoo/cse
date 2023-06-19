subroutine etkf ( x, m, n, delta, h, r, xobs, xda )

!*****************************************************************************80
!
!! EnKF (PO) data assimilation
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2022
!
!  Author:
!
!    Mao Ouyang
!
!  Input:
!
!    real ( kind = rk ) x(m, n): ensemble states.
!    integer m, n: size of ensemble matrix.
!    real ( kind = rk ) delta: inflation covariance.
!    real ( kind = rk ) h(n, n): observation matrix.
!    real ( kind = rk ) r(n, n): observation error covariance.
!    real ( kind = rk ) xobs(n): observation data.
!
!  Output:
!
!    real ( kind = rk ) xda(m, n): analysis ensemble states.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  integer, parameter :: lwmax = 1000

  integer n, m, i, info, ipiv(n), lwork

  real ( kind = rk ) x(m, n)
  real ( kind = rk ) xfm(n), xfp(m, n)
  real ( kind = rk ) yfm(n), yfp(m, n)
  real ( kind = rk ) delta
  real ( kind = rk ) h(n, n), r(n, n), rp(n, n)
  real ( kind = rk ) cm(m, n)
  real ( kind = rk ) pa(m, m), im(m, m), w(m)
  real ( kind = rk ) d(m, m), u(m, m)
  real ( kind = rk ) wap(m, m), wam(m)
  real ( kind = rk ) xobs(n), yf(n, m)
  real ( kind = rk ) xda(m, n), xam(n), xap(n,m)
  real ( kind = rk ) workp(lwmax)

  ! calculate the forecast ensemble mean and perturbation
  write (*, *) 'Forecast ensembles'
  write(*, '(5(f15.8))') x
  do i = 1, n
    xfm(i) = SUM(x(:, i)) / m
    xfp(:, i) = x(:, i) - xfm(i)
  end do
  write (*, *) 'Forecast ensemble mean'
  write (*, '(3(f15.8))') xfm

  ! add the inflation variance to the perturbation matrix
  xfp(:, :) = xfp(:, :) * (1.0 + delta)
  write (*, *) 'Forecast ensemble perturbation'
  write(*, '(5(f15.8))') xfp

  ! calculate observation ensemble mean and perturbation
  yf = MATMUL(h, TRANSPOSE(x))
  do i = 1, n
    yfm(i) = SUM(yf(i, :)) / m
    yfp(:, i) = yf(i, :) - yfm(i)
  end do
  write (*, *) 'Observation ensemble mean'
  write (*, '(3(f15.8))') yfm
  write (*, *) 'Observation ensemble perturbation'
  write (*, '(5(f15.8))') yfp

  ! calculate the kalman gain
  ! (1) calculate analysis error covariance
  !!! inverse of R 
  rp = r
  forall(i=1:n) rp(i, i) = 1.0 / rp(i, i)
  cm = MATMUL(yfp, rp)
  write (*, *) 'C matrix'
  write (*, '(5(f15.8))') cm
  im = 0
  forall(i=1:m) im(i, i) = 1
  pa = (m-1.0) * im + MATMUL(cm, TRANSPOSE(yfp))
  !!! SVD of Pa by LAPACK
  lwork = -1
  call DSYEV ('Vectors', 'Upper', m, pa, m, w, workp, lwork, info)
  lwork = min(lwmax, int(workp(1)))
  call DSYEV ('Vectors', 'Upper', m, pa, m, w, workp, lwork, info)
  write (*, *) 'SVD'
  write (*, '(5(f15.8))') w
  write (*, '(5(f15.8))') pa
  d = 0
  forall(i=1:m) d(i, i) = 1.0 / w(i)
  u = pa
  pa = MATMUL(MATMUL(u, d), TRANSPOSE(u)) 
  write (*, *) 'Analysis error covariance'
  write (*, '(5(f15.8))') pa
  ! (2) calculate the weight for deviation update
  forall(i=1:m) d(i, i) = 1.0 / SQRT(w(i))
  wap = SQRT(m - 1.0) * MATMUL(MATMUL(u, d), TRANSPOSE(u)) 
  write (*, *) 'Analysis ensemble deviation weight'
  write (*, '(5(f15.8))') wap
  ! (3) calculate the weight for ensemble mean updates
  wam = MATMUL(MATMUL(pa, cm), (xobs - yfm))
  write (*, *) 'Analysis ensemble mean weight'
  write (*, '(5(f15.8))') wam

  ! update the analysis ensemble
  xam = xfm + MATMUL(TRANSPOSE(xfp), wam)
  write (*, *) 'Analysis ensemble mean update'
  write (*, '(3(f15.8))') xam
  xap = MATMUL(TRANSPOSE(xfp), wap)
  write (*, *) 'Analysis ensemble deviation update'
  write (*, '(3(f15.8))') xap
  do i = 1 , n
    xda(:, i) = xam(i) + xap(i, :)
  end do
  write (*, *) 'Analysis update'
  write (*, '(5(f15.8))') xda

  return
end
