module mod_retardo
use precision_90
use constantes
use aceleraciones_lw_f90
use variables
use aceleraciones_lw_orden2_f90
implicit none

!==================== NOTACIÓN ====================
!LAS CANTIDADES PARA A SE EVALUAN ASI "(rA, rB, vA, vB, mA, mB)"
!LAS CANTIDADES PARA B SE EVALUAN ASI "(rB, rA, vB, vA, mA, mB)"

  !================== Parámetros iteración de Δ y τ′ ==================
real(kind=wp), parameter :: tol_dt       = 1.0e-12_wp
integer,      parameter :: maxit_dt      = 200
real(kind=wp), parameter :: tol_tau      = 1.0e-11_wp
integer,      parameter :: max_iter_tau  = 200

  ! --------- salidas persistentes del módulo ---------
  ! Orden 0 vistos por A y por B
real(kind=wp) :: dtB0_A=0.0_wp, rB0_A(3)=0.0_wp, vB0_A(3)=0.0_wp, aB0_A(3)=0.0_wp, aA0_A(3)=0.0_wp
real(kind=wp) :: dtA0_B=0.0_wp, rA0_B(3)=0.0_wp, vA0_B(3)=0.0_wp, aA0_B(3)=0.0_wp, aB0_B(3)=0.0_wp

  ! Orden 1
real(kind=wp) :: dtB1_A=0.0_wp, rB1_A(3)=0.0_wp, vB1_A(3)=0.0_wp, aB1_A(3)=0.0_wp, aA1_A(3)=0.0_wp
real(kind=wp) :: dtA1_B=0.0_wp, rA1_B(3)=0.0_wp, vA1_B(3)=0.0_wp, aA1_B(3)=0.0_wp, aB1_B(3)=0.0_wp

  ! Orden 2 (RK2 + τ′) —
real(kind=wp) :: dtp_A=0.0_wp, rB2_A(3)=0.0_wp, vB2_A(3)=0.0_wp, aA_stage_A(3)=0.0_wp, aB_stage_A(3)=0.0_wp, aA2_A(3)=0.0_wp
real(kind=wp) :: dtp_B=0.0_wp, rA2_B(3)=0.0_wp, vA2_B(3)=0.0_wp, aB_stage_B(3)=0.0_wp, aA_stage_B(3)=0.0_wp, aB2_B(3)=0.0_wp
real(kind=wp) :: rAtp_A(3)=0.0_wp, vAtp_A(3)=0.0_wp, rBtp_A(3)=0.0_wp, vBtp_A(3)=0.0_wp
real(kind=wp) :: rAtp_B(3)=0.0_wp, vAtp_B(3)=0.0_wp, rBtp_B(3)=0.0_wp, vBtp_B(3)=0.0_wp

real(kind=wp), save, public :: rA2_at_snap(3)=0._wp, vA2_at_snap(3)=0._wp
real(kind=wp), save, public :: rB2_at_snap(3)=0._wp, vB2_at_snap(3)=0._wp
real(kind=wp), save, public :: aA2_A_snap(3)=0._wp, aB2_B_snap(3)=0._wp
logical,      save, public :: snap_done_A = .false., snap_done_B = .false.

real(kind=wp)        :: tau0A, tau0B
contains
 ! ================== Δ0_A (ec. 13) ==================
function Delta0_A(rA, rB, vA)
real(kind=wp) :: Delta0_A
real(kind=wp) :: rA(3), rB(3), vA(3)
real(kind=wp) :: rAB(3), v2, fac, r2, rv_c2, disc

rAB   = rA - rB
v2    = dot3(vA, vA)
r2    = dot3(rAB, rAB)
fac   = 1.0_wp - v2/(cef*cef)
rv_c2 = dot3(rAB, vA)/(cef*cef)
disc  = fac*r2/(cef*cef) + rv_c2*rv_c2

Delta0_A = ( sqrt(max(0.0_wp, disc)) + rv_c2 ) / max(1.0e-30_wp, fac)
end function Delta0_A


! ================== Δ1_A (ec. 20) ==================

function Delta1_A(rA, rB, vA, vB, mA, mB)
real(kind=wp) :: Delta1_A
real(kind=wp) :: rA(3), rB(3), vA(3), vB(3), mA, mB
real(kind=wp) :: fA0(3), aA0(3)
real(kind=wp) :: dt, dt_new, rABc(3), v2, fac, r2, rv_c2, disc
integer :: it

  ! aA0 desde fuerza O(0)
fA0 = Force0_A(rA, vA, rB, vB, mA, mB)
aA0 = fA0 / mA


dt     = Delta0_A(rA, rB, vA)
dt_new = dt
do it = 1, maxit_dt
    rABc  = (rB - rA) - 0.5_wp*(dt*dt)*aA0
    v2    = dot3(vA, vA)
    r2    = dot3(rABc, rABc)
    fac   = 1.0_wp - v2/(cef*cef)
    rv_c2 = dot3(rABc, vA)/(cef*cef)
    disc  = fac*r2/(cef*cef) + rv_c2*rv_c2
    dt_new = ( sqrt(max(0.0_wp, disc)) + rv_c2 ) / max(1.0e-30_wp, fac)
    if (abs(dt_new - dt) <= tol_dt*max(1.0_wp, abs(dt))) exit
    dt = dt_new
end do
Delta1_A = dt_new
end function Delta1_A


! ================== Δp_A (ec. 57) ==================
function DeltaP_A(rA, rB, vA, vB, mA, mB)
real(kind=wp) :: DeltaP_A
real(kind=wp) :: rA(3), rB(3), vA(3), vB(3), mA, mB
real(kind=wp), parameter :: chi = 2.0_wp/3.0_wp
real(kind=wp), parameter :: alpha2 = 0.5_wp
real(kind=wp) :: dtB1, fA1(3), aA1(3)
real(kind=wp) :: fB0(3), aB0(3), rB_2o3(3)
real(kind=wp) :: dtp, dtp_new, rfn(3), v2, fac, r2, rv_c2, disc
integer :: it


fB0 = Force0_B(rB, vB, rA, vA, mA, mB)
aB0 = fB0 / mB

  ! Δt_B1 (herramienta)
dtB1 = Delta1_B(rB, rA, vB, vA, mA, mB)

  ! r_B(t - 2/3 Δt_B1) (ec. 55)
rB_2o3 = rB - (chi*dtB1)*vB + 0.5_wp*((chi*dtB1)**2)*aB0

  ! a_A1 (fuerza de orden 1)
fA1 = Force1_A(rA, vA, rB, vB, mA, mB)
aA1 = fA1 / mA


dtp     = chi * dtB1
dtp_new = dtp
do it = 1, maxit_dt
    rfn   = (rB_2o3 - rA) - alpha2*(dtp*dtp)*aA1   ! ec. (58)
    v2    = dot3(vA, vA)
    fac   = 1.0_wp - v2/(cef*cef)
    r2    = dot3(rfn, rfn)
    rv_c2 = dot3(rfn, vA)/(cef*cef)


    disc   = (rv_c2 + chi*dtB1)**2 + fac * ( (r2/(cef*cef)) - (chi*dtB1)**2 )
    dtp_new = ( rv_c2 + chi*dtB1 + sqrt(max(0.0_wp, disc)) ) / max(1.0e-30_wp, fac)

     if (abs(dtp_new - dtp) <= tol_dt*max(1.0_wp, abs(dtp))) exit
    dtp = dtp_new
end do
DeltaP_A = dtp_new
end function DeltaP_A



! ================== Δ0_B (ec. 13) ==================
function Delta0_B(rB, rA, vB)
real(kind=wp) :: Delta0_B
real(kind=wp) :: rB(3), rA(3), vB(3)
real(kind=wp) :: rBA(3), v2, fac, r2, rv_c2, disc

rBA   = rB - rA
v2    = dot3(vB, vB)
r2    = dot3(rBA, rBA)
fac   = 1.0_wp - v2/(cef*cef)
rv_c2 = dot3(rBA, vB)/(cef*cef)
disc  = fac*r2/(cef*cef) + rv_c2*rv_c2

Delta0_B = ( sqrt(max(0.0_wp, disc)) + rv_c2 ) / max(1.0e-30_wp, fac)
end function Delta0_B


! ================== Δ1_B (ec. 20) ==================
function Delta1_B(rB, rA, vB, vA, mA, mB)
real(kind=wp) :: Delta1_B
real(kind=wp) :: rB(3), rA(3), vB(3), vA(3), mA, mB
real(kind=wp) :: fB0(3), aB0(3)
real(kind=wp) :: dt, dt_new, rBAc(3), v2, fac, r2, rv_c2, disc
integer :: it

  ! aB0 desde fuerza O(0)
fB0 = Force0_B(rB, vB, rA, vA, mA, mB)
aB0 = fB0 / mB


dt     = Delta0_B(rB, rA, vB)
dt_new = dt
do it = 1, maxit_dt
    rBAc  = (rA - rB) - 0.5_wp*(dt*dt)*aB0
    v2    = dot3(vB, vB)
    r2    = dot3(rBAc, rBAc)
    fac   = 1.0_wp - v2/(cef*cef)
    rv_c2 = dot3(rBAc, vB)/(cef*cef)
    disc  = fac*r2/(cef*cef) + rv_c2*rv_c2
    dt_new = ( sqrt(max(0.0_wp, disc)) + rv_c2 ) / max(1.0e-30_wp, fac)
     if (abs(dt_new - dt) <= tol_dt*max(1.0_wp, abs(dt))) exit
    dt = dt_new
end do
Delta1_B = dt_new
end function Delta1_B


! ================== Δp_B (ec. 57) ==================
function DeltaP_B(rB, rA, vB, vA, mA, mB)
real(kind=wp) :: DeltaP_B
real(kind=wp) :: rB(3), rA(3), vB(3), vA(3), mA, mB
real(kind=wp), parameter :: chi = 2.0_wp/3.0_wp
real(kind=wp), parameter :: alpha2 = 0.5_wp
real(kind=wp) :: dtA1, fB1(3), aB1(3)
real(kind=wp) :: fA0(3), aA0(3), rA_2o3(3)
real(kind=wp) :: dtp, dtp_new, rfn(3), v2, fac, r2, rv_c2, disc
integer :: it


fA0 = Force0_A(rA, vA, rB, vB, mA, mB)
aA0 = fA0 / mA

  ! Δt_A1 (herramienta)
dtA1 = Delta1_A(rA, rB, vA, vB, mA, mB)

  ! r_A(t - 2/3 Δt_A1)
rA_2o3 = rA - (chi*dtA1)*vA + 0.5_wp*((chi*dtA1)**2)*aA0

  ! a_B1
fB1 = Force1_B(rB, vB, rA, vA, mA, mB)
aB1 = fB1 / mB


dtp     = chi * dtA1
dtp_new = dtp
do it = 1, maxit_dt
    rfn   = (rA_2o3 - rB) - alpha2*(dtp*dtp)*aB1
    v2    = dot3(vB, vB)
    fac   = 1.0_wp - v2/(cef*cef)
    r2    = dot3(rfn, rfn)
    rv_c2 = dot3(rfn, vB)/(cef*cef)

    disc   = (rv_c2 + chi*dtA1)**2 + fac * ( (r2/(cef*cef)) - (chi*dtA1)**2 )
    dtp_new = ( rv_c2 + chi*dtA1 + sqrt(max(0.0_wp, disc)) ) / max(1.0e-30_wp, fac)

     if (abs(dtp_new - dtp) <= tol_dt*max(1.0_wp, abs(dtp))) exit
    dtp = dtp_new
end do
DeltaP_B = dtp_new
end function DeltaP_B


  !==================== Fuerza O(0) SOBRE A ====================
function Force0_A(rA, vA, rB, vB, mA, mB) result(fA0)
real(kind=wp) :: rA(3), vA(3), rB(3), vB(3), mA, mB
real(kind=wp) :: fA0(3)
real(kind=wp) :: dtB0, rB0(3), vB0(3), aB0(3), aA0(3)

dtB0 = Delta0_B(rB, rA, vB)
rB0  = rB - dtB0*vB
vB0  = vB
aB0  = 0.0_wp
tau0A = tau0_LL(q1, mA)

  ! A <- B (obs = A, src = B)
call accel_LL_step_AB( rA, vA,                 &
                         rB0, vB0, aB0, q2, &
                         q1, mA, tau0A,  &
                         aA0 )
!call lw_accel_on_A(rA, vA, rB0, vB0, aB0, mA, aA0)
fA0 = mA * aA0
end function Force0_A

  !==================== Fuerza O(0) SOBRE B ====================
function Force0_B(rB, vB, rA, vA, mA, mB) result(fB0)
real(kind=wp) :: rB(3), vB(3), rA(3), vA(3), mA, mB
real(kind=wp) :: fB0(3)
real(kind=wp) :: dtA0, rA0(3), vA0(3), aA0(3), aB0(3)

dtA0 = Delta0_A(rA, rB, vA)
rA0  = rA - dtA0*vA
vA0  = vA
aA0  = 0.0_wp

tau0B = tau0_LL(q2, mB)

call accel_LL_step_AB( rB, vB,                 &
                        rA0, vA0,aA0, q1, &
                         q2, mB, tau0B,  &
                         aB0 )
!call lw_accel_on_B(rB, vB, rA0, vA0, aA0, mB, aB0)
fB0 = mB * aB0
end function Force0_B

  !==================== Fuerza O(1) SOBRE A ====================
function Force1_A(rA, vA, rB, vB, mA, mB) result(fA1)
real(kind=wp) :: rA(3), vA(3), rB(3), vB(3), mA, mB
real(kind=wp) :: fA1(3)
real(kind=wp) :: fB0(3), dtB1, aB0(3), rB1(3), vB1(3), aB1(3), aA1(3)

fB0 = Force0_B(rB, vB, rA, vA, mA, mB)
aB0 = fB0 / mB

dtB1 = Delta1_B(rB, rA, vB, vA, mA, mB)
rB1  = rB - dtB1*vB + 0.5_wp*(dtB1*dtB1)*aB0
vB1  = vB - dtB1*aB0
aB1  = aB0
tau0A = tau0_LL(q1, mA)
call accel_LL_step_AB( rA, vA,                 &
                         rB1, vB1, aB1, q2, &
                         q1, mA, tau0A,  &
                         aA1 )
!call lw_accel_on_A(rA, vA, rB1, vB1, aB1, mA, aA1)
fA1 = mA * aA1
end function Force1_A

  !==================== Fuerza O(1) SOBRE B ====================
function Force1_B(rB, vB, rA, vA, mA, mB) result(fB1)
real(kind=wp) :: rB(3), vB(3), rA(3), vA(3), mA, mB
real(kind=wp) :: fB1(3)
real(kind=wp) :: fA0(3), dtA1, aA0(3), rA1(3), vA1(3), aA1(3), aB1(3)

fA0 = Force0_A(rA, vA, rB, vB, mA, mB)
aA0 = fA0 / mA

dtA1 = Delta1_A(rA, rB, vA, vB, mA, mB)
rA1  = rA - dtA1*vA + 0.5_wp*(dtA1*dtA1)*aA0
vA1  = vA - dtA1*aA0
aA1  = aA0
tau0B = tau0_LL(q2, mB)
call accel_LL_step_AB( rB, vB,                 &
                         rA1, vA1, aA1, q1, &
                         q2, mB, tau0B,  &
                         aB1 )
!call lw_accel_on_B(rB, vB, rA1, vA1, aA1, mB, aB1)
fB1 = mB * aB1
end function Force1_B

  !==================== Fuerza O(2) SOBRE A (RK2 + τ′) ====================
function Force2_A(rA, vA, rB, vB, mA, mB) result(fA2)
real(kind=wp) :: rA(3), vA(3), rB(3), vB(3), mA, mB
real(kind=wp) :: fA2(3)
real(kind=wp), parameter :: chi = 2.0_wp/3.0_wp

real(kind=wp) :: fA0(3), fB0(3), fA1(3), fB1(3)
real(kind=wp) :: dtA1, dtB1, dtAp, dtBp, tB_star
real(kind=wp) :: aA0(3), aB0(3), aA1(3), aB1(3)
real(kind=wp) :: aA_stage(3), aB_stage(3), aA2(3)
real(kind=wp) :: rAtp(3), vAtp(3), rBtp(3), vBtp(3)
real(kind=wp) :: rB2_at(3), vB2_at(3)

fA0 = Force0_A(rA, vA, rB, vB, mA, mB);  aA0 = fA0/mA
fB0 = Force0_B(rB, vB, rA, vA, mA, mB);  aB0 = fB0/mB

dtA1 = Delta1_A(rA, rB, vA, vB, mA, mB)
dtB1 = Delta1_B(rB, rA, vB, vA, mA, mB)

fA1 = Force1_A(rA, vA, rB, vB, mA, mB);  aA1 = fA1/mA
fB1 = Force1_B(rB, vB, rA, vA, mA, mB);  aB1 = fB1/mB

dtAp = DeltaP_A(rA, rB, vA, vB, mA, mB)
dtBp = DeltaP_B(rB, rA, vB, vA, mA, mB)


 !----- lado A
 !if (dtA1 > 0.0_wp .and. t_now > 2*dtAp .and. dtB1 > 0.0_wp .and. t_now > 2*dtBp) then
  !  call set_dt1_A(t_now, dtA1, dtAp)   ! guarda teval_A = t_now - dtA1

!----- lado B

    !call set_dt1_B(t_now, dtB1, dtBp)   ! guarda teval_B = t_now - dtB1
 ! endif



rBtp = rB - dtBp*vB + 0.5_wp*dtBp*dtBp*aB1
vBtp = vB - dtBp*aB1
rAtp = rA - dtAp*vA + 0.5_wp*dtAp*dtAp*aA1
vAtp = vA - dtAp*aA1
! A <- B (obs = A, src = B)
tau0A = tau0_LL(q1, mA)
call accel_LL_step_AB(rA - chi*dtA1*vA, vA - chi*dtA1*aA1,                 &
                         rBtp, vBtp, aB1, q2, &
                         q1, mA, tau0A,  &
                         aA_stage )


call accel_LL_step_AB( rB - chi*dtB1*vB, vB - chi*dtB1*aB1,                 &
                         rAtp, vAtp, aA1, q1, &
                        q2, mB, tau0B,  &
                         aB_stage )
!call lw_accel_on_A( rA - chi*dtA1*vA,  &
                    !vA - chi*dtA1*aA1, &
                    !rBtp, vBtp, aB1, mA, aA_stage )

!call lw_accel_on_B( rB - chi*dtB1*vB,  &
                  !  vB - chi*dtB1*aB1, &
                   ! rAtp, vAtp, aA1, mB, aB_stage )

rB2_at = rB - dtB1*vB + 0.5_wp*dtB1*dtB1*aB1
vB2_at = vB - 0.25_wp*dtB1*aB1 - 0.75_wp*dtB1*aB_stage

call accel_LL_step_AB(rA, vA,                 &
                         rB2_at, vB2_at, aB1, q2, &
                         q1, mA, tau0A,  &
                         aA2 )
!call lw_accel_on_A( rA, vA, rB2_at, vB2_at, aB1, mA, aA2 )

fA2 = mA * aA2
aA2_A = aA2

!if (.not. snap_done_B) then
   ! rB2_at_snap = rB2_at
  !  vB2_at_snap = vB2_at
   ! aA2_A_snap  = aA2
   ! snap_done_B = .true.
!end if


end function Force2_A

  !==================== Fuerza O(2) SOBRE B (RK2 + τ′) ====================
function Force2_B(rB, vB, rA, vA, mA, mB) result(fB2)
real(kind=wp) :: rB(3), vB(3), rA(3), vA(3), mA, mB
real(kind=wp) :: fB2(3)
real(kind=wp), parameter :: chi = 2.0_wp/3.0_wp

real(kind=wp) :: fA0(3), fB0(3), fA1(3), fB1(3)
real(kind=wp) :: dtA1, dtB1, dtAp, dtBp, tA_star
real(kind=wp) :: aA0(3), aB0(3), aA1(3), aB1(3)
real(kind=wp) :: aA_stage(3), aB_stage(3), aB2(3)
real(kind=wp) :: rAtp(3), vAtp(3), rBtp(3), vBtp(3)
real(kind=wp) :: rA2_at(3), vA2_at(3)

fB0 = Force0_B(rB, vB, rA, vA, mA, mB);  aB0 = fB0/mB
fA0 = Force0_A(rA, vA, rB, vB, mA, mB);  aA0 = fA0/mA

dtB1 = Delta1_B(rB, rA, vB, vA, mA, mB)
dtA1 = Delta1_A(rA, rB, vA, vB, mA, mB)

fA1 = Force1_A(rA, vA, rB, vB, mA, mB);  aA1 = fA1/mA
fB1 = Force1_B(rB, vB, rA, vA, mA, mB);  aB1 = fB1/mB


dtAp = DeltaP_A(rA, rB, vA, vB, mA, mB)
dtBp = DeltaP_B(rB, rA, vB, vA, mA, mB)


rBtp = rB - dtBp*vB + 0.5_wp*dtBp*dtBp*aB1
vBtp = vB - dtBp*aB1
rAtp = rA - dtAp*vA + 0.5_wp*dtAp*dtAp*aA1
vAtp = vA - dtAp*aA1

tau0B = tau0_LL(q2, mB)
call accel_LL_step_AB( rB - chi*dtB1*vB, vB - chi*dtB1*aB1,                 &
                        rAtp, vAtp, aA1, q1, &
                         q2, mB, tau0B,  &
                         aB_stage )

call accel_LL_step_AB( rA - chi*dtA1*vA, vA - chi*dtA1*aA1,                 &
                        rBtp, vBtp, aB1, q2, &
                        q1, mA, tau0A,  &
                        aA_stage )
!call lw_accel_on_B( rB - chi*dtB1*vB,  &
!                      vB - chi*dtB1*aB1, &
!                      rAtp, vAtp, aA1, mB, aB_stage )

!call lw_accel_on_A( rA - chi*dtA1*vA,  &
 !                     vA - chi*dtA1*aA1, &
!                      rBtp, vBtp, aB1, mA, aA_stage )

rA2_at = rA - dtA1*vA + 0.5_wp*dtA1*dtA1*aA1
vA2_at = vA - 0.25_wp*dtA1*aA1 - 0.75_wp*dtA1*aA_stage


call accel_LL_step_AB( rB, vB,                 &
                         rA2_at, vA2_at, aA1, q1, &
                         q2, mB, tau0B,  &
                         aB2 )
!call lw_accel_on_B( rB, vB, rA2_at, vA2_at, aA1, mB, aB2 )

fB2 = mB * aB2
aB2_B = aB2
!if (.not. snap_done_A) then
    !rA2_at_snap = rA2_at
    !vA2_at_snap = vA2_at
   ! aB2_B_snap  = aB2
   ! snap_done_A = .true.
!end if
end function Force2_B

!===========================
!  ORDEN 0 – lado A
!===========================
subroutine retA0(rA, vA, rB, vB, mA, mB)
real(kind=wp) :: rA(3), vA(3), rB(3), vB(3), mA, mB
real(kind=wp) :: fA0(3)

! Δ0_B y estado fuente 0 para B (vista de A)
dtB0_A = Delta0_B(rB, rA, vB)
rB0_A  = rB - dtB0_A*vB
vB0_A  = vB
aB0_A  = 0.0_wp

! aA0 final vía fuerza O(0)
fA0    = Force0_A(rA, vA, rB, vB, mA, mB)
aA0_A  = fA0 / mA
end subroutine retA0

!===========================
!  ORDEN 0 – lado B
!===========================
subroutine retB0(rB, vB, rA, vA, mA, mB)
real(kind=wp) :: rB(3), vB(3), rA(3), vA(3), mA, mB
real(kind=wp) :: fB0(3)

! Δ0_A y estado fuente 0 para A (vista de B)
dtA0_B = Delta0_A(rA, rB, vA)
rA0_B  = rA - dtA0_B*vA
vA0_B  = vA
aA0_B  = 0.0_wp

  ! aB0 final vía fuerza O(0)
fB0    = Force0_B(rB, vB, rA, vA, mA, mB)
aB0_B  = fB0 / mB
end subroutine retB0

!===========================
!  ORDEN 1 – lado A
!===========================
subroutine retA1(rA, vA, rB, vB, mA, mB)
real(kind=wp) :: rA(3), vA(3), rB(3), vB(3), mA, mB
real(kind=wp) :: fB0_loc(3)

  ! Δ1_B “herramienta” (usa internamente Force0_B y Δ0)
dtB1_A = Delta1_B(rB, rA, vB, vA, mA, mB)


fB0_loc = Force0_B(rB, vB, rA, vA, mA, mB)
aB0_A   = fB0_loc / mB
rB1_A   = rB - dtB1_A*vB + 0.5_wp*(dtB1_A*dtB1_A)*aB0_A
vB1_A   = vB - dtB1_A*aB0_A
aB1_A   = aB0_A

  ! aA1 final vía fuerza O(1)
aA1_A = ( Force1_A(rA, vA, rB, vB, mA, mB) ) / mA
end subroutine retA1

!===========================
!  ORDEN 1 – lado B
!===========================
subroutine retB1(rB, vB, rA, vA, mA, mB)
real(kind=wp) :: rB(3), vB(3), rA(3), vA(3), mA, mB
real(kind=wp) :: fA0_loc(3)

  ! Δ1_A  (usa internamente Force0_A y Δ0)
dtA1_B = Delta1_A(rA, rB, vA, vB, mA, mB)

  !Para rA1/vA1, necesit0 aA0:
fA0_loc = Force0_A(rA, vA, rB, vB, mA, mB)
aA0_B   = fA0_loc / mA
rA1_B   = rA - dtA1_B*vA + 0.5_wp*(dtA1_B*dtA1_B)*aA0_B
vA1_B   = vA - dtA1_B*aA0_B
aA1_B   = aA0_B

  ! aB1 final vía fuerza O(1)
aB1_B = ( Force1_B(rB, vB, rA, vA, mA, mB) ) / mB
end subroutine retB1

!===========================
!  ORDEN 2 – lado A
!===========================
subroutine retA2(rA, vA, rB, vB, mA, mB)
real(kind=wp) :: rA(3), vA(3), rB(3), vB(3), mA, mB
dtB0_A = Delta0_B(rB, rA, vB)
dtA0_B = Delta0_A(rA, rB, vA)
  ! Δ1 de ambos (herramientas)
dtA1_B = Delta1_A(rA, rB, vA, vB, mA, mB)
dtB1_A = Delta1_B(rB, rA, vB, vA, mA, mB)

  ! Δp para A (usa Force1_A y Δ1_B internamente)
dtp_A  = DeltaP_A(rA, rB, vA, vB, mA, mB)


  ! aA2 final vía fuerza O(2) (usa internamente Δ1 y Δp)
aA2_A = ( Force2_A(rA, vA, rB, vB, mA, mB) ) / mA
end subroutine retA2

!===========================
!  ORDEN 2 – lado B
!===========================
subroutine retB2(rB, vB, rA, vA, mA, mB)
real(kind=wp) :: rB(3), vB(3), rA(3), vA(3), mA, mB

dtB0_A = Delta0_B(rB, rA, vB)
dtA0_B = Delta0_A(rA, rB, vA)
  ! Δ1 de ambos (herramientas)
dtB1_A = Delta1_B(rB, rA, vB, vA, mA, mB)
dtA1_B = Delta1_A(rA, rB, vA, vB, mA, mB)

  ! Δp para B (herramienta; usa Force1_B y Δ1_A internamente)
dtp_B  = DeltaP_B(rB, rA, vB, vA, mA, mB)


  ! aB2 final vía fuerza O(2) (usa internamente Δ1 y Δp)
aB2_B = ( Force2_B(rB, vB, rA, vA, mA, mB) ) / mB
end subroutine retB2


end module mod_retardo
