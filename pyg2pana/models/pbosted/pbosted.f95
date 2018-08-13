subroutine cal_xs_array(Z, A, E, Ep, theta, xs, N)
    integer*8, intent(in) :: Z, A
    double precision, dimension(N), intent(in) :: E, Ep, theta
    double precision, dimension(N), intent(out) :: xs
    integer :: N
    !f2py intent(hide) :: N=shape(Z,0)
    !f2py depend(N) xs
    integer :: i
    double precision :: xs_temp

    do i = 1, N
        call cal_xs_scalar(Z, A, E(i), Ep(i), theta(i), xs_temp)
        xs(i) = xs_temp
    end do

end subroutine cal_xs_array

subroutine cal_xs_scalar(Z, A, E, Ep, theta, xs)
    integer*8, intent(in) :: Z, A
    double precision, intent(in) :: E, Ep, theta
    double precision, intent(out) :: xs
    double precision :: Z1, A1
    double precision :: q2, w2, nu
    double precision :: ALPHA, M
    double precision :: sin2_theta_2, cos2_theta_2, tan2_theta_2
    double precision :: F1, F2, r, xs1, xs2, mott

    M = 0.93828  ! 0.93828 is used in F1F209.f
    ALPHA = 1 / 137.0388
    Z1 = Z
    A1 = A

    nu = E - Ep
    
    sin2_theta_2 = sin(abs(theta) / 2.0)**2
    cos2_theta_2 = 1.0 - sin2_theta_2
    tan2_theta_2 = sin2_theta_2 / cos2_theta_2

    q2 = 4.0 * E * Ep * sin2_theta_2
    w2 = M**2 + 2.0 * M * nu - q2

    mott = (ALPHA / (2 * E * sin2_theta_2))**2 * cos2_theta_2 * 389.379

    call F1F2IN09(Z1, A1, q2, w2, F1, F2, r)
    xs1 = mott * (2.0 / M * F1 * tan2_theta_2 + F2 / nu)

    call F1F2QE09(Z1, A1, q2, w2, F1, F2)
    xs2 = mott * (2.0 / M * F1 * tan2_theta_2 + F2 / nu)

    xs = (xs1 + xs2) / 1000.0  !ub/MeV-sr

end subroutine cal_xs_scalar
