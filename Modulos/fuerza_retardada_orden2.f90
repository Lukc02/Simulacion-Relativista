module aceleraciones_lw_orden2_f90
use precision_90
use constantes
use aceleraciones_lw_f90
implicit none


contains

!------------------------------------------------------------
function tau0_LL(q, m)
! Escala de reacción radiativa (Landau–Lifshitz):
!   tau0 = (2 * k1_ef * q^2) / (3 * m * c^3)
real(kind=wp), intent(in) :: q, m
real(kind=wp) :: tau0
real(kind=wp) :: c3
real(kind=wp) :: tau0_LL
c3   = cef*cef*cef
tau0_LL = (2.0_wp * k1_ef *q*q) / (3.0_wp * m * c3)
end function tau0_LL
  !------------------------------------------------------------

  !------------------------------------------------------------
function gamma_from_v(v)
  ! γ = 1 / sqrt(1 - |v|^2/c^2)
real(kind=wp), intent(in) :: v(3)
real(kind=wp) :: gamma
real(kind=wp) :: beta2
real(kind=wp) :: gamma_from_v


beta2 = dot3(v,v) / cef**2
gamma_from_v = 1.0_wp / sqrt(1.0_wp - beta2)
end function gamma_from_v
  !------------------------------------------------------------

  !------------------------------------------------------------
function accel_from_fourforce(v, f, m)
  ! A.33: m a = (1/γ^2) f − (1/γ^2) (f·v) v / c^2
real(kind=wp), intent(in) :: v(3), f(3), m
real(kind=wp) :: gamma, invg2, fv, c2
real(kind=wp), dimension(3) :: accel_from_fourforce

c2    = cef*cef
gamma = gamma_from_v(v)
invg2 = 1.0_wp/(gamma*gamma)
fv    = dot3(f, v)

accel_from_fourforce = (invg2/m) * ( f - (fv/c2)*v )
end function accel_from_fourforce
 !------------------------------------------------------------
subroutine fL_spatial(q, v, E, B, fL)
! f^i = q F^{i\nu} u_\nu  →  f_vec = q*γ*(E + v×B)
real(kind=wp), intent(in)  :: q, v(3), E(3), B(3)
real(kind=wp), intent(out) :: fL(3)
real(kind=wp) :: gamma, vxB(3)
gamma = gamma_from_v(v)
call cross_product(v, B, vxB)
fL = q*gamma*( E + vxB )
end subroutine fL_spatial

!------------------------------------------------------------
subroutine lw_fields_velonly_tdot(x_obs, v_obs, r_tr, v_tr, a_tr, q, dEdt, dBdt)
! Derivada temporal en t del término "velocidad" de LW:
!   E_vel = k q (1 - beta^2) (n - beta) / (kappa^3 R^2)
! pero escrita aquí como  (γ^{-2})(n - beta)/(kappa^3 R^2).
! *** IMPORTANTE: todas las magnitudes con subíndice "tr" están evaluadas en t_r(t),
! y por cadena temporal: d/dt = (1/kappa) d/dt_r sobre ellas. ***

real(kind=wp), intent(in)  :: x_obs(3), v_obs(3), r_tr(3), v_tr(3), a_tr(3), q
real(kind=wp), intent(out) :: dEdt(3), dBdt(3)

real(kind=wp) :: Rv(3), n(3), R, invR2, invR3
real(kind=wp) :: beta(3), betad_t(3), beta2, gamma, inv_g2
real(kind=wp) :: nb, kap, kap2, kap3, kap4
real(kind=wp) :: RdotVec(3), Rdot, ndot(3), ndot_beta, g2dot_t, scale_dot
real(kind=wp) :: E0(3), nxE(3)

! --- Geometría observador-fuente retardada ---
Rv = x_obs - r_tr
R  = sqrt( dot3(Rv,Rv) )
n      = Rv / R
invR2  = 1.0_wp/(R*R)
invR3  = invR2 / R

! --- Estado fuente (retardado) ---
beta   =  v_tr / cef
beta2  =  dot3(beta,beta)
gamma  =  1.0_wp / sqrt( 1.0_wp - beta2)
inv_g2 =  1.0_wp / (gamma*gamma)

nb        = dot3(n,beta)
kap       = 1.0_wp - nb
kap2      = kap*kap
kap3      = kap2*kap
kap4      = kap3*kap

! --- Reglas de cadena temporal (FIX(κ)):
! d r_tr / dt = (1/kap) v_tr,  d v_tr / dt = (1/kap) a_tr
! => betad_t = d beta / dt = (1/c) d v_tr / dt = (1/kap) a_tr / c
betad_t = (1.0_wp/kap) * (a_tr/cef)

! dR/dt = n·(v_obs - d r_tr/dt) = n·(v_obs - v_tr/kap)
RdotVec   = v_obs - (1.0_wp/kap)*v_tr      ! FIX(κ)
Rdot      = dot3(n, RdotVec)

! dn/dt = (RdotVec - n*(n·RdotVec)) / R
ndot      = (RdotVec - n*Rdot) / R

! d/dt (n·beta) = (dn/dt)·beta + n·(dbeta/dt)
ndot_beta = dot3(ndot, beta) + dot3(n, betad_t)

! d/dt (γ^{-2}) = d/dt (1 - beta^2) = -2 beta·betad_t   (FIX(κ) ya incluido en betad_t)
g2dot_t   = -2.0_wp * dot3(beta, betad_t)

! --- Campo "velocidad" base para dB/dt ---
E0 = k1_ef*q * (n - beta) * (invR2/kap3) * inv_g2

! --- d/dt [ (γ^{-2}) (kappa^{-3}) (R^{-2}) ] (regla del producto)
scale_dot = g2dot_t*(invR2/kap3) + &
            inv_g2*( 3.0_wp*ndot_beta*(invR2/kap4) - 2.0_wp*Rdot*(invR3/kap3) )

! --- dE_vel/dt = k q [ d/dt(n - beta)*(...) + (n - beta)*d/dt((...)) ]
dEdt = k1_ef*q * ( (ndot - betad_t)*(invR2/kap3)*inv_g2 + (n - beta)*scale_dot )

! --- dB_vel/dt (SI):  B = (1/c) n×E  =>  dB/dt = (ndot×E0 + n×dE/dt)/c
call cross_product(ndot, E0, nxE)
call cross_product(n,    dEdt, dBdt)
dBdt = (nxE + dBdt) / cef
end subroutine lw_fields_velonly_tdot



subroutine fRR_LL_spatial(tau0, q, m, v, E, B, dEdt, dBdt, fRR)
! Landau–Lifshitz (2º orden) según ec. (68):
! f_RR = tau0 * c^3 * { q γ (Ė + v×Ḃ)  + (q/m)[ f0L E + fL×B ]
!                       - (γ/(m c^2)) ( f0L^2 - |fL|^2 ) v }
! donde: fL = q γ (E + v×B),  f0L = (q γ / c) (E·v)

real(kind=wp), intent(in)  :: tau0, q, m
real(kind=wp), intent(in)  :: v(3), E(3), B(3), dEdt(3), dBdt(3)
real(kind=wp), intent(out) :: fRR(3)

real(kind=wp) :: gamma, vxdB(3), tmp(3)
real(kind=wp) :: fL(3), f0L, Fdot_u(3), F_times_f(3), F_proj(3)
real(kind=wp) :: fL2, c2

c2    = cef*cef
gamma = gamma_from_v(v)

! (i) f_L (espacial) y f_L^0 (temporal)
call fL_spatial(q, v, E, B, fL)
f0L = (q*gamma/cef) * dot3(E, v)                ! f^0_L = γ q (E·v)/c

! (ii) q γ (Ė + v×Ḃ)
call cross_product(v, dBdt, vxdB)
Fdot_u = q*gamma * ( dEdt + vxdB )

! (iii) (q/m)[ f0L E + fL×B ]
call cross_product(fL, B, tmp)
F_times_f = (q/m) * ( f0L*E + tmp )

! (iv) Término cuadrático/proyector:  − (γ/(m c^2)) ( f0L^2 − |fL|^2 ) v
fL2 = dot3(fL, fL)
F_proj = - (gamma/(m*c2)) * ( -f0L*f0L + fL2 ) * v

! (v) Ensamble final
fRR =  tau0 * ( Fdot_u + F_times_f + F_proj )
end subroutine fRR_LL_spatial


subroutine accel_LL_step_AB(x_obs, v_obs, r_src_tr, v_src_tr, a_src_tr, q_src, q_obs, m_obs, tau0_obs, a_obs)
real(kind=wp), intent(in)  :: x_obs(3), v_obs(3)
real(kind=wp), intent(in)  :: r_src_tr(3), v_src_tr(3), a_src_tr(3)
real(kind=wp), intent(in)  :: q_src, q_obs, m_obs, tau0_obs
real(kind=wp), intent(out) :: a_obs(3)

real(kind=wp) :: E(3), B(3), dEdt(3), dBdt(3), fL(3), fRR(3), fTot(3)

! 1) Campos LW COMPLETOS
call lw_fields(x_obs, r_src_tr, v_src_tr, a_src_tr, q_src, E, B)

! 2) “Punto” SOLO sobre el término Coulomb
call lw_fields_velonly_tdot(x_obs, v_obs, r_src_tr, v_src_tr, a_src_tr, q_src, dEdt, dBdt)

! 3) Fuerza de Lorentz espacial
call fL_spatial(q_obs, v_obs, E, B, fL)

! 4) Corrección LL (orden 2)
call fRR_LL_spatial(tau0_obs, q_obs, m_obs, v_obs, E, B, dEdt, dBdt, fRR)

! 5) Total y A.33
fTot = fL + fRR
a_obs = accel_from_fourforce(v_obs, fTot, m_obs)
end subroutine accel_LL_step_AB

end module aceleraciones_lw_orden2_f90
