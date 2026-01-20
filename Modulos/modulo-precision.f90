module precision_90
implicit none

integer, parameter :: s_precision = selected_real_kind(p=precision(0.0))
integer, parameter :: d_precision = selected_real_kind(p=precision(0.0_s_precision)+1)
integer, parameter :: tq_precision = selected_real_kind(p=precision(0.0_d_precision)+1)
integer, parameter :: q_precision = selected_real_kind(p=precision(0.0_tq_precision)+1)

!integer, parameter :: wp = s_precision
!integer, parameter :: wp = d_precision
integer, parameter :: wp = tq_precision
!integer, parameter :: wp = q_precision

end module precision_90
