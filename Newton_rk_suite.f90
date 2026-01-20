PROGRAM glob_dyn_time
use precision_90
use rksuite_90
use constantes
use variables
use aceleraciones_lw_f90
use derivada_retardo
use modulo_lagran4
use funciones
implicit none

integer :: k
real(kind=wp) ::  t_cpu_start, t_cpu_end, t_cpu
logical :: consistency_done
real(kind=wp) :: dtA1, dtB1, dtAp, dtBp             ! retardos
real(kind=wp) :: t_eval_A, t_eval_B ,t_cumplido                ! tiempos de evaluación
real(kind=wp) :: rA0(3), vA0(3), rB0(3), vB0(3)
real(kind=wp) :: t4_num(4)
real(kind=wp) :: rA4_num(3,4), vA4_num(3,4)
real(kind=wp) :: rB4_num(3,4), vB4_num(3,4)
real(kind=wp) :: rA_interp_num(3), vA_interp_num(3)
real(kind=wp) :: rB_interp_num(3), vB_interp_num(3)
real(kind=wp) :: rel_dR_A(3), rel_dV_A(3), rel_dA_A(3)
real(kind=wp) :: rel_dR_B(3), rel_dV_B(3), rel_dA_B(3)
real(kind=wp) :: denom_rA, denom_vA, denom_aA
real(kind=wp) :: denom_rB, denom_vB, denom_aB
real(kind=wp) :: dR_A_num(3), dV_A_num(3), dR_B_num(3), dV_B_num(3)
real(kind=wp) :: rA(3), vA(3), rB(3), vB(3)
real(kind=wp) :: dA_A_num(3), dA_B_num(3)
real(kind=wp) :: aA1_now(3), aB1_now(3), aB_model_interp(3), aA_model_interp(3)
real(kind=wp) :: dtA1_test, dtB1_test, dtAp_test, dtBp_test
real(kind=wp) :: aA1(3), aB1(3)
integer :: totf, step_cost, num_succ
real(kind=wp) :: true_error, waste, max_error, t_max_error, dummy ,mcheps ,dwarf , t_want, tiempoini, tiempofin, th_tol
real(kind=wp) :: unit_len , unit_mass , unit_time, unit_chr
real(kind=wp) :: t, h, En0, En0_eV, En, Lx, Ly, Lz, Px, Py, Pz
real(kind=wp) :: Lx_0, Ly_0, Lz_0
real(kind=wp) :: Px_0, Py_0, Pz_0
real(kind=wp) :: X_CM_0, Y_CM_0, Z_CM_0
real(kind=wp) :: t_prev
real(kind=wp) :: R_CM_0(3)
real(kind=wp) :: a, ecc, v0, r0,p
real(kind=wp) ::  dy(12)
real(kind=wp) :: X_CM, Y_CM, Z_CM
real(kind=wp) :: R_CM(3)
real(kind=wp) :: unit_energy_SI, En_eV
real(kind=wp) :: dtA0_ini, dtB0_ini, dtA1_ini, dtB1_ini, dtAp_ini, dtBp_ini
real(kind=wp) :: accel_ini_A(3), accel_ini_B(3)
character(len=50) :: hache, hache0,metodo
character(len=90) :: archivo
logical :: error_assess0
type(rk_comm_real_1d) :: comm
integer :: flag


namelist /parametrosinic/ t_start , t_end
namelist /constantfis/ km, GN, msol, cluz, epsilon0, mu0, hplnck, hbarex, eelect, melect, unosalf, jz0
namelist /parametroscmp/ tol, metodo, orden
!namelist /paramstresh_Newton/ tresh

!open  (unit = 11, file = 'datosent/paramstresh_Newton.in', status = 'old' )
!read  (unit = 11, *) (tresh(i) , i=1,12)
open (unit= 10, file = 'datosent/paramscmp_Newton.in', status = 'old' )
read (unit= 10, nml = parametroscmp)

open(unit = 22, file='datosent/parametrosinic-N.in', status='old')
read(unit = 22, nml=parametrosinic )
open(unit = 24,  file='datosent/constantesfis.in', status='old')
read(unit = 24, nml=constantfis )

close(unit = 10)
close(unit = 11)
close(unit = 22)
close(unit = 24)

open(unit = 101 , file='datossal/salidas-New.f' , status='unknown')

pi = 4.0_wp*atan(1.0_wp)

!write (hache0 ,'( es10 .3) ') h

!hache = 'h = ' // trim ( adjustl ( hache0 ) )

write(* ,*) ' '
write(101 ,*) ' '

write(* ,*) ' metodo =', metodo
write(101,*) ' metodo =', metodo
write(*,*) 'tol = ',tol
write(101,*) 'tol = ', tol
write(*,*) ' orden = ', orden
write(101,*) ' orden = ', orden
write(*,*) ' precicion =', wp
write(101,*) ' precision =', wp
write(*,*) ' cluz =', cluz
write(101,*) 'cluz = ', cluz
write(*,*) ' GN =', GN
write(101,*) 'GN = ', GN
write(*,*) ' msol =', msol
write(101,*) 'msol = ', msol
write(*,*) ' cluz =', cluz
write(101,*) 'cluz = ', cluz

write(*,*) 'hplnck =', hplnck
write(*,*) 'hbarex =', hbarex
write(*,*) 'eelect =', eelect
write(*,*) 'melect =', melect
write(*,*) 'unosalf =', unosalf
elprmtt = (eelect**2)*unosalf/(2._wp*hplnck*cluz)
write(*,*) 'elprmtt =', elprmtt
write(*,*) 'k1_elect =', k1_elect   !! Jackson

k1_elex = 1._wp/(4._wp*pi*elprmtt)  !! Definicion exacta
write(*,*) 'k1_elex =', k1_elex
write(*,*) '10e-7 c**2 =', (cluz**2)*(1.e-7_wp)
abhr = unosalf*(hplnck)/(2._wp*pi*melect*cluz)
write(*,*) 'abhr =', abhr

write(101,*) 'hplnck =', hplnck
write(101,*) 'hbarex =', hbarex
write(101,*) 'eelect =', eelect
write(101,*) 'melect =', melect
write(101,*) 'unosalf =', unosalf
write(101,*) 'elprmtt =', elprmtt
write(101,*) 'k1_elect =', k1_elect   !! Jackson
write(101,*) 'k1_elex =', k1_elex
write(101,*) '10e-7 c**2 =', (cluz**2)*(1.e-7_wp)
write(101,*) 'abhr =', abhr

! Eleccion de unidades
unit_len = abhr !Radio de Bohr
unit_mass = melect !Masa del electron
unit_time = unit_len/cluz !Tiempo de vuelo para un radio de Bohr
unit_chr = eelect !Carga del electron

cef = cluz*(unit_time/unit_len)
Gef = GN * (unit_mass/unit_len**3)  !Constante de la gravedad en km**3/(s**2 * msol)
Gef_c2 = GN * (unit_mass/((cluz**2)*(unit_len**3))) !Constante de la gravedad en km

write(*,*) ' unit_len =', unit_len
write(101,*) 'unit_len = ', unit_len
write(*,*) ' unit_mass =', unit_mass
write(101,*) 'unit_mass = ', unit_mass
write(*,*) ' unit_time =', unit_time
write(101,*) 'unit_time = ', unit_time
write(*,*) ' unit_chr =', unit_chr
write(101,*) 'unit_chr = ', unit_chr

k1_ef = k1_elex * unit_chr**2 * unit_time**2/(unit_len**3 * unit_mass)

jz0_ef = jz0/unosalf

write(*,*) 'cef =', cef
write(101,*) 'cef =', cef
write(*,*) 'k1_ef =', k1_ef
write(101,*) 'k1_ef =', k1_ef
write(*,*) ' Gef (km**3/(s**2 * msol)) =' , Gef
write(101,*)'Gef (km**3/(s**2 * msol)) = ' , Gef
write(*,*) 'jz0_ef =', jz0_ef
write(101,*) 'jz0_ef =', jz0_ef

q1 =  1._wp
q2 = -1._wp

M1 = 1836.152673426_wp
M2 = 1._wp


if (jz0_ef == 0) then
    ecc = 0.4_wp
    a = 7.5_wp !En unidades del radio de Bohr

    r0 = a * (1.0_wp - ecc)

    periodo = 2._wp*pi*sqrt((M1*M2)/(k1_ef*abs(q1*q2)*(M1+M2))) * sqrt(a**3)
else
    ecc = 0.0_wp
    a   = jz0_ef**2 * (M1 + M2) /(M1*M2*k1_ef*abs(q1*q2)*(1-ecc**2))

    r0 = a * (1.0_wp - ecc)

    periodo = 2._wp*pi*sqrt((M1*M2)/(k1_ef*abs(q1*q2)*(M1+M2))) * sqrt(a**3)
end if

y_start(1:6) = (/ (r0* M2)/(M1+M2)  , 0.0_wp, 0.0_wp, &
        -(r0 * M1)/(M1+M2) , 0.0_wp, 0.0_wp /)

y_start(7:12) = &
         (/ 0.0_wp, sqrt( k1_ef*abs(q1*q2) *M2* (1.0_wp + ecc) / (a * (1.0_wp - ecc)*M1*(M1+M2) )), 0.0_wp, &
         0.0_wp, -sqrt( (k1_ef*abs(q1*q2) *M1* (1.0_wp + ecc)) / (a * (1.0_wp - ecc)*M2*(M1+M2) )) , 0.0_wp /)

dtA0_ini = Delta0_A(y_start(1:3), y_start(4:6), y_start(7:9))
dtA1_ini = Delta1_A(y_start(1:3), y_start(4:6), y_start(7:9), y_start(10:12), M1, M2)
dtAp_ini = DeltaP_A(y_start(1:3), y_start(4:6), y_start(7:9), y_start(10:12), M1, M2)
dtB0_ini = Delta0_B(y_start(4:6), y_start(1:3), y_start(10:12))
dtB1_ini = Delta1_B(y_start(4:6), y_start(1:3), y_start(10:12), y_start(7:9), M1, M2)
dtBp_ini = DeltaP_B(y_start(4:6), y_start(1:3), y_start(10:12), y_start(7:9), M1, M2)

accel_ini_A = Force2_A(y_start(1:3), y_start(7:9), y_start(4:6), y_start(10:12), M1, M2)/M1
accel_ini_B = Force2_B(y_start(4:6), y_start(10:12), y_start(1:3), y_start(7:9), M1, M2)/M2

write(*,*) ' dtA0_ini =' , dtA0_ini
write(101,*) 'dtA0_ini = ' , dtA0_ini

write(*,*) ' dtA1_ini =' , dtA1_ini
write(101,*) 'dtA1_ini = ' , dtA1_ini

write(*,*) ' dtAp_ini =' , dtAp_ini
write(101,*) 'dtAp_ini = ' , dtAp_ini

write(*,*) ' dtB0_ini =' , dtB0_ini
write(101,*) 'dtB0_ini = ' , dtB0_ini

write(*,*) ' dtB1_ini =' , dtB1_ini
write(101,*) 'dtB1_ini = ' , dtB1_ini

write(*,*) ' dtBp_ini =' , dtBp_ini
write(101,*) 'dtBp_ini = ' , dtBp_ini

write(*,*) ' accel_ini_A =' , accel_ini_A
write(101,*) 'accel_ini_A = ' , accel_ini_A

write(*,*) ' accel_ini_B =' , accel_ini_B
write(101,*) 'accel_ini_B = ' , accel_ini_B

write(*,*) ' ecc =' , ecc
write(101,*) 'ecc = ' , ecc

write(*,*) ' a =' , a
write(101,*) 'a  = ' , a

write(*,*) ' r0  =' , r0
write(101,*) 'r0  = ' , r0

!write (* ,*) ' v0  =' , v0
!write (101 ,*) 'v0  = ' , v0

write(*,*) ' M1 (msol) =' , M1
write(101,*) 'M1 (msol) = ' , M1

write(*,*) ' M2 (msol) =' , M2
write(101,*) 'M2 (msol) = ' , M2

write(*,*) ' q1  =' , q1
write(101,*) 'q1  = ' , q1

write(*,*) ' q2  =' , q2
write(101,*) 'q2  = ' , q2

write(*,*) ' Gef_c2 =' , Gef_c2
write(101,*) 'Gef_c2 = ' , Gef_c2

write(*,*) ' periodo =' , periodo
write(101,*) 'periodo = ' , periodo

write(*,*) ' t_end/periodo =' , t_end/(2._wp*pi*sqrt((a**3)*(M1*M2)/((k1_ef*abs(q1*q2)*(M1+M2)))))
write(101,*) 't_end/periodo = ' , t_end/(2._wp*pi*sqrt((a**3)*(M1*M2)/((k1_ef*abs(q1*q2)*(M1+M2)))))

write(*,*) ' t_start = ' , t_start
write(*,*) ' t_end = ' , t_end
write(*,*) 'y_start = ' , y_start
write(101,*) ' t_start = ' , t_start
write(101,*) ' t_end = ' , t_end

write(101,*) ' y_start = ' , y_start

write(*,*) ' tol =',tol

mcheps = EPSILON (dummy) !It is the smallest value of x such that 1.+x is not equal to 1.
write(*,*)'      mcheps = ',mcheps
dwarf = TINY (dummy)
write(*,*)'       dwarf = ',dwarf !smallest positive number of a numeric data type that can be
                               !represented without loss of precision

write(*,*)   ' sqrt(dwarf) = ',SQRT(dwarf)
write(101,*) ' mcheps = ',mcheps ,'It is the smallest value of x such that 1.+x is not equal to 1.'
write(101,*) ' dwarf = ', dwarf , 'smallest positive number of a numeric data type that can be'
write(101,*) '                               represented without loss of precision'

write(101,*) 'sqrt(dwarf) = ',SQRT(dwarf)

write(* ,*) ' '
write(101 ,*) ' '

write(*,*) 'derivt(0._wp, y_start) = ' , derivt (0._wp, y_start)
write(10,*) 'derivt(0._wp, y_start) = ' , derivt (0._wp, y_start)

call constantes_movimiento(y_start, En0, Lx_0, Ly_0, Lz_0, Px_0, Py_0, Pz_0, X_CM_0, Y_CM_0, Z_CM_0)

R_CM_0 = (/ X_CM_0, Y_CM_0, Z_CM_0 /)

thres(1:6) = (/ 1e-20_wp, 1e-20_wp, 1e-20_wp, 1e-20_wp, 1e-20_wp, 1e-20_wp /)
thres(7:12) = (/ 1e-20_wp, 1e-20_wp, 1e-20_wp, 1e-20_wp, 1e-20_wp, 1e-20_wp /)
error_assess0 = .false.

t_prev = t_start

write(* ,*) ' En0 = ' , En0
write(101 ,*) ' En0 = ' , En0
write(*,*)  'Lx_0 = ' , Lx_0
write(101,*) 'Lx_0 = ' , Lx_0
write(*,*)  'Ly_0 = ' , Ly_0
write(101,*)'Ly_0 = ' , Ly_0
write(*,*)  'Lz_0 = ' , Lz_0
write(101,*)'Lz_0 = ' , Lz_0
write(*,*)  'Px_0 = ' , Px_0
write(101,*)'Px_0 = ' , Px_0
write(*,*)  'Py_0 = ' , Py_0
write(101,*)'Py_0 = ' , Py_0
write(*,*)  'Pz_0 = ' , Pz_0
write(101,*)'Pz_0 = ' , Pz_0
write(*,*)  'R_CM_0 = ', R_CM_0
write(101,*) 'R_CM_0 = ', R_CM_0
write(*,*)   'th_tol =',th_tol
write(101,*) 'th_tol =',th_tol

write(*,*) 'Evaluaciones_derivada = ', n_derivt_calls
write(101,*) 'Evaluaciones_derivada = ', n_derivt_calls

archivo = 'datossal/EM_posiciones_ret.dat'
open(unit =10 , file = archivo)

archivo = 'datossal/EM_velocidades_ret.dat'
open(unit =11, file = archivo)

archivo = 'datossal/EM_aceleraciones_ret.dat'
open(unit =12, file = archivo)

archivo = 'datossal/EM_MA_ret.dat'
open(unit =30, file = archivo)

archivo = 'datossal/EM_momento_ret.dat'
open(unit =40, file = archivo)

archivo = 'datossal/EM_ene_ret.dat'
open(unit =50, file = archivo)

archivo = 'datossal/EM_CM_ret.dat'
open(unit =60, file = archivo)

write(10,*) 't', ' X_1x ' , ' X_1y ' , ' X_1z ' , ' X_2x ' , ' X_2y ' , ' X_2z '

write(11,*) 't',  ' V_1x ', ' V_1y ' , ' V_1z ' , ' V_2x ' ,' V_2y ' , ' V_2z '

write(12,*) 't',  ' a_1x ', ' a_1y ' , ' a_1z ' , ' a_2x ' ,'a_2y ' , ' a_2z'

write(30,*) 't' , 'Lx', 'Ly', 'Lz', '(Lz_0 - Lz)/Lz_0'

write(40,*) 't' , 'Px', 'Py', 'Pz', '(Py - Py_0)/Py_0'

write(50,*) 't' , 'En', ' (En - En0)/En0 '

write(60,*) 't', 'X_CM', 'Y_CM', 'Z_CM' , '(X_CM_0 - X_CM)/X_CM_0', '(Y_CM_0 - Y_CM)/Y_CM_0', 'Z_CM_0 - Z_CM)/Z_CM_0'


write(10,'(E29.8,6(E29.13))') t_start, y_start(1), y_start(2), y_start(3), y_start(4), y_start(5), y_start(6)
write(11,'(E29.8,6(E29.13))') t_start, y_start(7), y_start(8), y_start(9),y_start(10), y_start(11), y_start(12)
!write(12,'(E29.8,6(E29.13))') t_now, dy(7), dy(8), dy(9),&
                                    !dy(10), dy(11), dy(12)

write(30,'(E29.8,4(E29.13))') t_start , Lx_0 , Ly_0, Lz_0, 0._wp
write(40,'(E29.8,4(E29.13))') t_start , Px_0 , Py_0, Pz_0, 0._wp
write(50,'(E29.8,2(E29.13))') t_start , En0,  0._wp
write(60,'(E29.8,6(E29.13))') t_start, X_CM_0, Y_CM_0, Z_CM_0 , 0._wp , 0._wp, 0._wp


call setup(comm,t_start,y_start,t_end,tol,thres,method=metodo,task='step', &
           error_assess=error_assess0)
!           error_assess=.false.)
!           error_assess=.true.)

error_assess0 = .true.
! ==> INICIO MEDICIÓN TIEMPO
call cpu_time(t_cpu_start)
consistency_done = .false.

do
    call step_integrate(comm,derivt,t_now,y_now=y_now,flag=flag)

    call constantes_movimiento(y_now, En, Lx, Ly, Lz, Px, Py, Pz, X_CM, Y_CM, Z_CM)
    dy = derivt(t_now,y_now)

    write(10,'(E29.8,6(E29.13))') t_now, y_now(1), y_now(2), y_now(3), y_now(4), y_now(5), y_now(6)
    write(11,'(E29.8,6(E29.13))') t_now, y_now(7), y_now(8), y_now(9),y_now(10), y_now(11), y_now(12)
    write(12,'(E29.8,6(E29.13))') t_now, dy(7), dy(8), dy(9),&
                                    dy(10), dy(11), dy(12)

    write(30,'(E29.8,4(E29.13))') t_now , Lx , Ly, Lz, ABS((Lz_0 - Lz)/Lz_0)
    write(40,'(E29.8,4(E29.13))') t_now , Px , Py, Pz,  ABS((Py - Py_0)/Py_0)
    write(50,'(E29.8,2(E29.13))') t_now , En,   ABS((-En + En0)/En0)
    write(60,'(E29.8,6(E29.13))') t_now, X_CM, Y_CM, Z_CM , (-X_CM_0  +X_CM), (-Y_CM_0 + Y_CM), &
                                        (-Z_CM_0 + Z_CM)

    dtA1_test = Delta1_A(y_now(1:3), y_now(4:6), y_now(7:9), y_now(10:12), M1, M2)
    dtB1_test = Delta1_B(y_now(4:6), y_now(1:3), y_now(10:12), y_now(7:9), M1, M2)
    dtAp_test = DeltaP_A(y_now(1:3), y_now(4:6), y_now(7:9), y_now(10:12), M1, M2)
    dtBp_test = DeltaP_B(y_now(4:6), y_now(1:3), y_now(10:12), y_now(7:9), M1, M2)

        ! ================== CARGA INICIAL DE 4 PUNTOS (k = 1..4) ==================
    if (k >= 1 .and. k <= 4) then
        write(*,'(I6,1X,ES29.8,9(1X,ES29.13))') k, t_now, y_now(1:3), y_now(7:9)
        write(*,'(I6,1X,ES29.8,9(1X,ES29.13))') k, t_now, y_now(4:6), y_now(10:12)

        t4_num(    k-0) = t_now
        rA4_num(:, k-0) = y_now(1:3);  vA4_num(:, k-0) = y_now(7:9)
        rB4_num(:, k-0) = y_now(4:6);  vB4_num(:, k-0) = y_now(10:12)
    end if


    if (.not. consistency_done) then
  ! --- 1) Retardos de 1er y 2do orden en el presente ---
        dtA1 = Delta1_A(y_now(1:3), y_now(4:6), y_now(7:9), y_now(10:12), M1, M2)
        dtB1 = Delta1_B(y_now(4:6), y_now(1:3), y_now(10:12), y_now(7:9), M1, M2)
        dtAp = DeltaP_A(y_now(1:3), y_now(4:6), y_now(7:9), y_now(10:12), M1, M2)
        dtBp = DeltaP_B(y_now(4:6), y_now(1:3), y_now(10:12), y_now(7:9), M1, M2)

        if (t_now > 5.0_wp*dtAp .and. t_now > 5.0_wp*dtBp) then
    ! --- 2) Congelar PRESENTE del paso que cumple la condición ---
            rA0 = y_now(1:3);  vA0 = y_now(7:9)
            rB0 = y_now(4:6);  vB0 = y_now(10:12)

    ! --- 3) Tiempos de evaluación (para la trayectoria numérica) ---
            t_eval_A = t_now - dtA1
            t_eval_B = t_now - dtB1
            snap_done_A = .false.
            snap_done_B = .false.

            write(*,'(/,A)') '=============================================='
            write(*,'(A)')   '    Condición de tiempo retardado cumplida'
            write(*,'(A,I16)') '       k     = ', k
            write(*,'(A,E29.8)') '       t_now     = ', t_now
            write(*,'(A,E29.8)') '       deltat_A     = ', dtAp_test
            write(*,'(A,E29.8)') '       deltat_B     = ', dtBp_test
            write(*,'(A,E29.8)') '       deltat_A1     = ', dtA1_test
            write(*,'(A,E29.8)') '       deltat_B1     = ', dtB1_test
            write(*,'(A,E29.8)') '       t_eval_A  = ', t_eval_A
            write(*,'(A,E29.8)') '       t_eval_B  = ', t_eval_B
            write(*,'(A)')   '    Se habilita el test de consistencia global'
            write(*,'(A)')   '=============================================='

            write(101,'(/,A)') '=============================================='
            write(101,'(A)')   '   Condición de tiempo retardado cumplida'
            write(101,'(A,I16)') '       k     = ', k
            write(101,'(A,E29.8)') '       t_now     = ', t_now
            write(101,'(A,E29.8)') '       delta_t_A     = ', dtAp_test
            write(101,'(A,E29.8)') '       delta_t_B     = ', dtBp_test
            write(101,'(A,E29.8)') '       delta_t_A1     = ', dtA1_test
            write(101,'(A,E29.8)') '       delta_t_B1     = ', dtB1_test
            write(101,'(A,E29.8)') '       t_eval_A  = ', t_eval_A
            write(101,'(A,E29.8)') '       t_eval_B  = ', t_eval_B
            write(101,'(A)')   '   Se habilita el test de consistencia global'
            write(101,'(A)')   '=============================================='

            call retA2(rA0, vA0, rB0, vB0, M1, M2)
            call retB2(rB0, vB0, rA0, vA0, M1, M2)

            !==============================================
            write(*,'(/,A)') '=============================================='
            write(*,'(A)')   '    Cantidades del modelo (lado A)'
            write(*,'(A,I6)') '       k          = ', k
            write(*,'(A,E29.8)') '       t_now      = ', t_now
            write(*,'(A,3(1X,E22.12))') '       rA2_at_snap = ', rA2_at_snap
            write(*,'(A,3(1X,E22.12))') '       vA2_at_snap = ', vA2_at_snap
            write(*,'(A,3(1X,E22.12))') '       aA2_A_snap  = ', aA2_A_snap
            write(*,'(A,E29.8)') '       t_eval_B   = ', t_eval_B
            write(*,'(A)')   '=============================================='

            !==============================================
            write(*,'(/,A)') '=============================================='
            write(*,'(A)')   '    Cantidades del modelo (lado B)'
            write(*,'(A,I6)') '       k          = ', k
            write(*,'(A,E29.8)') '       t_now      = ', t_now
            write(*,'(A,3(1X,E22.12))') '       rB2_at_snap = ', rB2_at_snap
            write(*,'(A,3(1X,E22.12))') '       vB2_at_snap = ', vB2_at_snap
            write(*,'(A,3(1X,E22.12))') '       aB2_B_snap  = ', aB2_B_snap
            write(*,'(A,E29.8)') '       t_eval_A   = ', t_eval_A
            write(*,'(A)')   '=============================================='

            !==============================================
            write(101,'(/,A)') '=============================================='
            write(101,'(A)')   '   Cantidades del modelo (lado A)'
            write(101,'(A,I6)') '       k          = ', k
            write(101,'(A,E29.8)') '       t_now      = ', t_now
            write(101,'(A,3(1X,E22.12))') '       rA2_at_snap = ', rA2_at_snap
            write(101,'(A,3(1X,E22.12))') '       vA2_at_snap = ', vA2_at_snap
            write(101,'(A,3(1X,E22.12))') '       aA2_A_snap  = ', aA2_A_snap
            write(101,'(A,E29.8)') '       t_eval_B   = ', t_eval_B
            write(101,'(A)')   '=============================================='

            !==============================================
            write(101,'(/,A)') '=============================================='
            write(101,'(A)')   '   Cantidades del modelo (lado B)'
            write(101,'(A,I6)') '       k          = ', k
            write(101,'(A,E29.8)') '       t_now      = ', t_now
            write(101,'(A,3(1X,E22.12))') '       rB2_at_snap = ', rB2_at_snap
            write(101,'(A,3(1X,E22.12))') '       vB2_at_snap = ', vB2_at_snap
            write(101,'(A,3(1X,E22.12))') '       aB2_B_snap  = ', aB2_B_snap
            write(101,'(A,E29.8)') '       t_eval_A   = ', t_eval_A
            write(101,'(A)')   '=============================================='
            consistency_done = .true.

        end if
    end if

    k = k + 1
   if (t_now==t_end .or. flag > 3) exit
end do
call cpu_time(t_cpu_end)
t_cpu  = t_cpu_end - t_cpu_start

write(*,*) 'Tiempo CPU =', t_cpu


! ==== Interpolación  ====
call lagrange4_vec3(t4_num, vB4_num, t_eval_B, vB_interp_num)
call lagrange4_vec3(t4_num, rB4_num, t_eval_B, rB_interp_num)
call lagrange4_vec3(t4_num, rA4_num, t_eval_A, rA_interp_num)
call lagrange4_vec3(t4_num, vA4_num, t_eval_A, vA_interp_num)

aA1_now = Force1_A(rA0, vA0, rB0, vB0, M1, M2) / M1
aB1_now = Force1_B(rB0, vB0, rA0, vA0, M1, M2) / M2

call lw_accel_on_B(rB0, vB0, rA_interp_num, vA_interp_num, aA1_now, M2, aB_model_interp)

call lw_accel_on_A(rA0, vA0, rB_interp_num, vB_interp_num, aB1_now, M1, aA_model_interp)

! === Comparación de consistencia para la partícula A ===
denom_rA = sqrt( rA_interp_num(1)**2 + rA_interp_num(2)**2 + rA_interp_num(3)**2 )
denom_vA = sqrt( vA_interp_num(1)**2 + vA_interp_num(2)**2 + vA_interp_num(3)**2 )
denom_aA = sqrt( aA2_A_snap(1)**2 + aA2_A_snap(2)**2 + aA2_A_snap(3)**2 )

denom_rB = sqrt( rB_interp_num(1)**2 + rB_interp_num(2)**2 + rB_interp_num(3)**2 )
denom_vB = sqrt( vB_interp_num(1)**2 + vB_interp_num(2)**2 + vB_interp_num(3)**2 )
denom_aB = sqrt( aB2_B_snap(1)**2 + aB2_B_snap(2)**2 + aB2_B_snap(3)**2 )

! ===== Diferencias ABSOLUTAS =====
dR_A_num = rA_interp_num - rA2_at_snap
dV_A_num = vA_interp_num - vA2_at_snap
dR_B_num = rB_interp_num - rB2_at_snap
dV_B_num = vB_interp_num - vB2_at_snap

! ===== Relativos componente a componente, dividiendo por el MÓDULO =====
rel_dR_A(1) = abs(dR_A_num(1)) / denom_rA
rel_dR_A(2) = abs(dR_A_num(2)) / denom_rA
rel_dR_A(3) = abs(dR_A_num(3)) / denom_rA

rel_dV_A(1) = abs(dV_A_num(1)) / denom_vA
rel_dV_A(2) = abs(dV_A_num(2)) / denom_vA
rel_dV_A(3) = abs(dV_A_num(3)) / denom_vA


rel_dR_B(1) = abs(dR_B_num(1)) / denom_rB
rel_dR_B(2) = abs(dR_B_num(2)) / denom_rB
rel_dR_B(3) = abs(dR_B_num(3)) / denom_rB

rel_dV_B(1) = abs(dV_B_num(1)) / denom_vB
rel_dV_B(2) = abs(dV_B_num(2)) / denom_vB
rel_dV_B(3) = abs(dV_B_num(3)) / denom_vB

! Diferencias de aceleración (modelo presente con fuente retardada vs modelo retardado):
dA_A_num = aA_model_interp - aA2_A_snap
dA_B_num = aB_model_interp - aB2_B_snap

! Denominadores: módulo (norma) del vector de referencia
!denom_aA = sqrt(aA2_A_snap(1)**2 + aA2_A_snap(2)**2 + aA2_A_snap(3)**2)
!denom_aB = sqrt(aB2_B_snap(1)**2 + aB2_B_snap(2)**2 + aB2_B_snap(3)**2)

! Relativos componente a componente (sin tolerancias)
rel_dA_A(1) = abs(dA_A_num(1)) / denom_aA
rel_dA_A(2) = abs(dA_A_num(2)) / denom_aA
rel_dA_A(3) = abs(dA_A_num(3)) / denom_aA

rel_dA_B(1) = abs(dA_B_num(1)) / denom_aB
rel_dA_B(2) = abs(dA_B_num(2)) / denom_aB
rel_dA_B(3) = abs(dA_B_num(3)) / denom_aB


! ==========
write(*,'(/,A)')                 '===== Consistencia (A) ====='
write(*,'(A,3(1X,ES22.12))')     'delta_r_A  =', dR_A_num
write(*,'(A,3(1X,ES22.12))')     'delta_v_A  =', dV_A_num
write(*,'(A,3(1X,ES22.12))')     'delta_a_A  =', dA_A_num
write(*,'(A,3(1X,ES22.12))')     'delta_r_A / |r_A_interp| =', rel_dR_A
write(*,'(A,3(1X,ES22.12))')     'delta_v_A / |v_A_interp| =', rel_dV_A
write(*,'(A,3(1X,ES22.12))')     'delta_a_A / |a_A_modelR| =', rel_dA_A   ! denom = |aA2_A_snap|

write(*,'(/,A)')                 '===== Consistencia (B) ====='
write(*,'(A,3(1X,ES22.12))')     'delta_r_B  =', dR_B_num
write(*,'(A,3(1X,ES22.12))')     'delta_v_B  =', dV_B_num
write(*,'(A,3(1X,ES22.12))')     'delta_a_B  =', dA_B_num
write(*,'(A,3(1X,ES22.12))')     'delta_r_B / |r_B_interp| =', rel_dR_B
write(*,'(A,3(1X,ES22.12))')     'delta_v_B / |v_B_interp| =', rel_dV_B
write(*,'(A,3(1X,ES22.12))')     'delta_a_B / |a_B_modelR| =', rel_dA_B   ! denom = |aB2_B_snap|
write(*,*) ' '
write(*,'(A,3(1X,ES22.12))') 'rA_interp_num = ', rA_interp_num
write(*,'(A,3(1X,ES22.12))') 'vA_interp_num = ', vA_interp_num
write(*,'(A,3(1X,ES22.12))') 'rB_interp_num = ', rB_interp_num
write(*,'(A,3(1X,ES22.12))') 'vB_interp_num = ', vB_interp_num
write(*,'(A,3(1X,ES22.12))') 'aB_interp_num = ', aB_model_interp
write(*,'(A,3(1X,ES22.12))') 'aA_interp_num = ', aA_model_interp

! ===== Prints compactos (mismo hilo en pantalla) =====
write(101,'(/,A)')                 '===== Consistencia (A) ====='
write(101,'(A,3(1X,ES22.12))')     'delta_r_A  =', dR_A_num
write(101,'(A,3(1X,ES22.12))')     'delta_v_A  =', dV_A_num
write(101,'(A,3(1X,ES22.12))')     'delta_a_A  =', dA_A_num
write(101,'(A,3(1X,ES22.12))')     'delta_r_A / |r_A_interp| =', rel_dR_A
write(101,'(A,3(1X,ES22.12))')     'delta_v_A / |v_A_interp| =', rel_dV_A
write(101,'(A,3(1X,ES22.12))')     'delta_a_A / |a_A_modelR| =', rel_dA_A   ! denom = |aA2_A_snap|

write(101,'(/,A)')                 '===== Consistencia (B) ====='
write(101,'(A,3(1X,ES22.12))')     'delta_r_B  =', dR_B_num
write(101,'(A,3(1X,ES22.12))')     'delta_v_B  =', dV_B_num
write(101,'(A,3(1X,ES22.12))')     'delta_a_B  =', dA_B_num
write(101,'(A,3(1X,ES22.12))')     'delta_r_B / |r_B_interp| =', rel_dR_B
write(101,'(A,3(1X,ES22.12))')     'delta_v_B / |v_B_interp| =', rel_dV_B
write(101,'(A,3(1X,ES22.12))')     'delta_a_B / |a_B_modelR| =', rel_dA_B   ! denom = |aB2_B_snap|
write(101,*) ' '
write(101,'(A,3(1X,ES22.12))') 'rA_interp_num = ', rA_interp_num
write(101,'(A,3(1X,ES22.12))') 'vA_interp_num = ', vA_interp_num
write(101,'(A,3(1X,ES22.12))') 'rB_interp_num = ', rB_interp_num
write(101,'(A,3(1X,ES22.12))') 'vB_interp_num = ', vB_interp_num
write(101,*) 'Tiempo de CPU', t_cpu

close(101)

close(10)
close(30)
close(40)
close(50)
close(60)


END PROGRAM glob_dyn_time
