module derivada_retardo
use precision_90
use constantes
use variables
use aceleraciones_lw_f90
use mod_retardo
use funciones
implicit none

contains


function derivt(t_now, y)
real(kind=wp), intent(in) :: t_now
real(kind=wp), intent(in) :: y(:)
real(kind=wp)             :: derivt(size(y))


real(kind=wp) :: rA(3), rB(3), vA(3), vB(3)

real(kind=wp) :: rA_ret(3), vA_ret(3), aA(3)
real(kind=wp) :: rB_ret(3), vB_ret(3), aB(3)


rA = y(1:3)
rB = y(4:6)
vA = y(7:9)
vB = y(10:12)

select case (orden)
  case (0)
    ! Orden 0: usa Δ0 y fuerza O(0) de ambos lados
    call retA0(rA, vA, rB, vB, M1, M2)
    call retB0(rB, vB, rA, vA, M1, M2)
    aA = aA0_A    ! publicados por mod_retardo
    aB = aB0_B

  case (1)
    ! Orden 1: usa Δ1 con f0 y fuerza O(1)
    call retA1(rA, vA, rB, vB, M1, M2)
    call retB1(rB, vB, rA, vA, M1, M2)
    aA = aA1_A
    aB = aB1_B

  case (2)
    ! Orden 2: usa Δp con f1 y fuerza O(2)
    call retA2(rA, vA, rB, vB, M1, M2)
    call retB2(rB, vB, rA, vA, M1, M2)
    aA = aA2_A
    aB = aB2_B
end select

!call retardos_orden2_sincronizado(rA, vA, rB, vB, M1, M2,  &
 !                                 rA_ret, vA_ret, aA_ret,   &
 !                                 rB_ret, vB_ret, aB_ret)


derivt(1:3)   = vA
derivt(4:6)   = vB
derivt(7:9)   = aA
derivt(10:12) = aB

n_derivt_calls = n_derivt_calls + 1

end function derivt

end module derivada_retardo
