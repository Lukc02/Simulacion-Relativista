module modulo_lagran4
use precision_90
implicit none
private
public :: lagrange4_scalar, lagrange4_vec3

contains

subroutine lagrange4_scalar(t, y, te, ye)
  ! Lagrange baricéntrico con 4 nodos (estable y simple)
real(kind=wp), intent(in)  :: t(4), y(4), te
real(kind=wp), intent(out) :: ye
real(kind=wp) :: w(4), denom, numer
integer  :: i, j

    ! Si te coincide con un nodo, retornamos y(i) directamente
  do i=1,4
    if (te == t(i)) then
      ye = y(i); return
    end if
  end do

    ! Pesos baricéntricos: w_i = 1 / prod_{j≠i} (t_i - t_j)
  do i=1,4
    w(i) = 1.0_wp
    do j=1,4
      if (j /= i) w(i) = w(i) / (t(i) - t(j))
    end do
  end do

    ! Evaluación: p(te) = (sum w_i y_i /(te - t_i)) / (sum w_i /(te - t_i))
numer = 0._wp; denom = 0._wp
  do i=1,4
    numer = numer + w(i)*y(i)/(te - t(i))
    denom = denom + w(i)/(te - t(i))
  end do

ye = numer / denom
end subroutine lagrange4_scalar

subroutine lagrange4_vec3(t4, Y, te, yout)
real(kind=wp), intent(in)  :: t4(4)         ! tiempos
real(kind=wp), intent(in)  :: Y(3,4)        ! valores (3 componentes × 4 puntos)
real(kind=wp), intent(in)  :: te
real(kind=wp), intent(out) :: yout(3)

integer  :: ic
real(kind=wp) :: ytmp(4)

  do ic = 1, 3
     ! recolección explícita de la FILA ic en un buffer contiguo
     ytmp(1) = Y(ic,1)
     ytmp(2) = Y(ic,2)
     ytmp(3) = Y(ic,3)
     ytmp(4) = Y(ic,4)

     call lagrange4_scalar(t4, ytmp, te, yout(ic))
  end do
end subroutine lagrange4_vec3
end module modulo_lagran4
