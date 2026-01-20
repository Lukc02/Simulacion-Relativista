module aceleraciones_lw_f90
use precision_90
use constantes
implicit none
contains

subroutine cross_product(a, b, c)
real(kind=wp), intent(in)  :: a(3), b(3)
real(kind=wp), intent(out) :: c(3)

c(1) = a(2)*b(3) - a(3)*b(2)
c(2) = a(3)*b(1) - a(1)*b(3)
c(3) = a(1)*b(2) - a(2)*b(1)

end subroutine cross_product

function dot3(a,b)
real(kind=wp), intent(in) :: a(3), b(3)
real(kind=wp)             :: dot3

dot3 = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

end function dot3

subroutine lw_fields(x_obs, r_tr, v_tr, a_tr, q, E, B)
! Campos de Liénard–Wiechert
real(kind=wp), intent(in)  :: x_obs(3), r_tr(3), v_tr(3), a_tr(3), q
real(kind=wp), intent(out) :: E(3), B(3)
real(kind=wp) :: Rv(3), n(3), beta(3), dbeta(3), tmp(3), rad(3)
real(kind=wp) :: R, beta2, gamma, nb, den, prefC, prefR

Rv = x_obs - r_tr
R  = sqrt( dot3(Rv,Rv) )
 if (R == 0.0_wp) then
 E = 0.0_wp; B = 0.0_wp
   return
 end if

n = Rv / R

beta  = v_tr / cef
dbeta = a_tr / cef
beta2 = dot3(beta,beta)
gamma = 1.0_wp / sqrt( max(0.0_wp, 1.0_wp - beta2) )
nb    = dot3(n, beta)
den   = 1.0_wp - nb

! E = E_coulomb + E_radiation
prefC = q * k1_ef / ( (gamma*gamma) * (den**3) * (R*R) )
E = prefC * ( n - beta )

call cross_product( (n - beta), dbeta, tmp )
call cross_product( n, tmp, rad )

prefR = q * k1_ef / ( cef * (den**3) * R )
E = E + prefR * rad

! B = n × E / c
call cross_product(n, E, B)

B = B / cef

end subroutine lw_fields

subroutine lw_accel_on_A(xA, vA, rB_ret, vB_ret, aB_ret, mA, aA)
real(kind=wp), intent(in)  :: xA(3), vA(3), rB_ret(3), vB_ret(3), aB_ret(3), mA
real(kind=wp), intent(out) :: aA(3)
real(kind=wp) :: E(3), B(3), vxB(3), v2, gammaA, dotEv, fac

call lw_fields(xA, rB_ret, vB_ret, aB_ret, q2, E, B)

v2     = dot3(vA,vA)
gammaA = 1.0_wp / sqrt( max(0.0_wp, 1.0_wp - v2/(cef*cef)) )
call cross_product(vA, B, vxB)
dotEv  = dot3(E, vA)

aA = (q1 / (mA * gammaA))*(E + vxB - ((dotEv/(cef*cef)) * vA))

end subroutine lw_accel_on_A

subroutine lw_accel_on_B(xB, vB, rA_ret, vA_ret, aA_ret, mB, aB)
real(kind=wp), intent(in)  :: xB(3), vB(3), rA_ret(3), vA_ret(3), aA_ret(3), mB
real(kind=wp), intent(out) :: aB(3)
real(kind=wp) :: E(3), B(3), vxB(3), v2, gammaB, dotEv, fac

call lw_fields(xB, rA_ret, vA_ret, aA_ret, q1, E, B)

v2     = dot3(vB,vB)
gammaB = 1.0_wp / sqrt( max(0.0_wp, 1.0_wp - v2/(cef*cef)) )
call cross_product(vB, B, vxB)
dotEv  = dot3(E, vB)

aB = (q2 / (mB * gammaB))*(E + vxB - ((dotEv/(cef*cef)) * vB))

end subroutine lw_accel_on_B

end module aceleraciones_lw_f90
