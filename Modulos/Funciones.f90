module funciones
use precision_90
use constantes
implicit none
contains

subroutine constantes_movimiento(q, En, Lx, Ly, Lz, Px, Py, Pz, X_CM, Y_CM, Z_CM)

real(kind=wp), dimension(1:12), intent(in) :: q

real(kind=wp), intent(out) :: En, Lx, Ly, Lz, Px, Py, Pz
real(kind=wp), intent(out) :: X_CM, Y_CM, Z_CM

real(kind=wp), dimension(3) :: r1, r2, v1, v2, n
real(kind=wp) :: u1_2, u2_2, gamma1, gamma2, r12, pot


r1 = q(1:3)
r2 = q(4:6)
v1 = q(7:9)
v2 = q(10:12)


u1_2 = dot_product(v1, v1)
u2_2 = dot_product(v2, v2)


gamma1 = 1._wp / sqrt(1._wp - u1_2 / cef**2)
gamma2 = 1._wp / sqrt(1._wp - u2_2 / cef**2)


r12 = sqrt(sum( (r1-r2)**2 ))
n   = (r1-r2)/r12

! --- Energía potencial
pot = k1_ef * q1 * q2 * ( 1._wp/r12 + &
     ((1._wp/r12)*dot_product(v1,v2) + 1._wp / (2._wp * r12**3 * cef**2) ) * &
     dot_product(v1, r1-r2) * dot_product(v2, r1-r2) )

! --- Energía total ---
En = (gamma1 - 1._wp) * M1 * cef**2 + &
     (gamma2 - 1._wp) * M2 * cef**2 + pot

! --- Momento angular total ---
Lx = gamma1*M1 * (r1(2) * v1(3) - r1(3) * v1(2)) + &
     gamma2*M2 * (r2(2) * v2(3) - r2(3) * v2(2))
Ly = gamma1*M1 * (r1(3) * v1(1) - r1(1) * v1(3)) + &
     gamma2*M2 * (r2(3) * v2(1) - r2(1) * v2(3))
Lz = gamma1*M1 * (r1(1) * v1(2) - r1(2) * v1(1)) + &
     gamma2*M2 * (r2(1) * v2(2) - r2(2) * v2(1))

! --- Momento lineal total ---
Px = gamma1*M1 * v1(1) + gamma2*M2 * v2(1)
Py = gamma1*M1 * v1(2) + gamma2*M2 * v2(2)
Pz = gamma1*M1 * v1(3) + gamma2*M2 * v2(3)

! --- Centro de masa ---
X_CM = (M1 * r1(1) + M2 * r2(1)) / (M1 + M2)
Y_CM = (M1 * r1(2) + M2 * r2(2)) / (M1 + M2)
Z_CM = (M1 * r1(3) + M2 * r2(3)) / (M1 + M2)
end subroutine constantes_movimiento

end module funciones
