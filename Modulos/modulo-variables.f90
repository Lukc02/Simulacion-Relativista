module variables
use precision_90

real(kind=wp) :: t_start,  t_end,  tol, t_now
integer       :: orden
integer       :: n_derivt_calls = 0
real(kind=wp), dimension(12) :: y_start, y_now, thres, y_maxvals,assessed_error


end module variables
