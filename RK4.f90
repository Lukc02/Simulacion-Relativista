MODULE rk4

use precsn
use derivada

contains


subroutine rk4_sys(t,w,h,w_h)

implicit none

real(kind=wp), intent(in)   :: t,h
real(kind=wp),dimension(:), intent(in)   :: w
real(kind=wp),dimension(:), intent(out)   :: w_h
real(kind=wp),dimension(size(w))    :: k1,k2,k3,k4

k1 = derivat(t,w)

!print*, "k1 =", k1

k2 = derivat(t + 0.5_wp*h, w + 0.5_wp*h*k1)

!print*, "k2 =", k2

k3 = derivat(t + 0.5_wp*h, w + 0.5_wp*h*k2)

!print*, "k3 =", k3

k4 = derivat(t + h, w + h*k3)

!print*, "k4 =", k4

w_h = w + h*(k1 + 2._wp*k2 + 2._wp*k3 + k4)/6._wp

return
end subroutine

END MODULE
