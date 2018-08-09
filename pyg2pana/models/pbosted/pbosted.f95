subroutine cal_xs(Z, A, E, Ep, theta, xs, N)
    integer*8, intent(in) :: Z, A
    double precision, dimension(N), intent(in) :: E, Ep, theta
    double precision, dimension(N), intent(out) :: xs
    integer :: N
    !f2py intent(hide) :: N=shape(Z,0)
    !f2py depend(N) xs
    integer :: i
    double precision :: Z1, A1
    double precision :: q2, w2, M, nu
    double precision :: sin_theta_2, cos_theta_2
    double precision :: F1, F2, r, xs1, xs2

    M = 0.93825
    Z1 = Z
    A1 = A

    do i = 1, N
        nu = E(i) - Ep(i)
        
        sin_theta_2 = sin(abs(theta(i)) / 2.0)
        cos_theta_2 = sqrt(1 - sin_theta_2**2)

        q2 = 4.0 * E(i) * Ep(i) * sin_theta_2**2
        w2 = M**2 + 2.0 * M * nu - q2

        call F1F2IN09(Z1, A1, q2, w2, F1, F2, r)

        xs1 = (2.0 / 137.0 * Ep(i) / q2 * cos_theta_2)**2  !mott
        xs1 = xs1 * (2.0 / M * F1 * (sin_theta_2 / cos_theta_2)**2 + F2 / nu)
        xs1 = xs1 * 389.379

        call F1F2QE09(Z1, A1, q2, w2, F1, F2)
        
        xs2 = (2.0 / 137.0 * Ep(i) / q2 * cos_theta_2)**2  !mott
        xs2 = xs2 * (2.0 / M * F1 * (sin_theta_2 / cos_theta_2)**2 + F2 / nu)
        xs2 = xs2 * 389.379

        xs(i) = (xs1 + xs2) / 1000.0  !ub/MeV-sr
    end do
end subroutine cal_xs
