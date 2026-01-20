PROGRAM newton
use precision_90
use constantes
use funciones
use derivada
use rksuite_90

implicit none
real(kind=wp) :: true_error, waste, max_error, t_max_error, dummy ,mcheps ,dwarf , t_want, tiempoini, tiempofin, th_tol
real(kind=wp) :: unit_len , unit_mass , unit_time, unit_chr
real(kind=wp) :: t, h, En0, En0_eV, En, Lx, Ly, Lz, Px, Py, Pz, En0_normalized, darwin_corr
real(kind=wp) :: Lx_0, Ly_0, Lz_0
real(kind=wp) :: Px_0, Py_0, Pz_0
real(kind=wp) :: X_CM_0, Y_CM_0, Z_CM_0
real(kind=wp) :: t_prev, emax_r, emax_v, rms_r, rms_v
integer :: ncount, ulog
real(kind=wp), dimension(3) :: R_CM_0, n_start
real(kind=wp) :: a, ecc, v0, r0, precession_deg,p
real(kind=wp),dimension(12) ::  dy
real(kind=wp) :: rA_num(3), vA_num(3), rB_num(3), vB_num(3)
real(kind=wp) :: aB_from_interp(3), aA_from_interp(3)
real(kind=wp) :: erA, evA, erB, evB, eaB, eaA
real(kind=wp) :: y_der(12), y_int(12)
integer :: i
real(kind=wp) :: X_CM, Y_CM, Z_CM
real(kind=wp), dimension(3) :: R_CM
real(kind=wp) :: unit_energy_SI, En_eV
real(kind=wp), dimension(3)  :: v01, v02
real(kind=wp) :: v01_sqr, v02_sqr, gamma_start1, gamma_start2, r12_start
real(kind=wp), dimension(3) :: r_int_A, r_int_B, v_int_A, v_int_B, a_int_A, a_int_B
real(kind=wp), parameter :: eps_rel = 1.0e-20_wp
real(kind=wp) :: drA(3), dvA(3), daB(3)
real(kind=wp) :: nrpA, nvpA, napB
real(kind=wp) :: drB(3), dvB(3), daA(3)
real(kind=wp) :: nrpB, nvpB, napA
real(kind=wp) :: dtA0_ini, dtB0_ini, dtA1_ini, dtB1_ini, dtAp_ini, dtBp_ini
real(kind=wp) :: accel_ini_A(3), accel_ini_B(3)
character(len=50) :: hache, hache0,metodo
character(len=90) :: archivo
logical :: error_assess0
type(rk_comm_real_1d) :: comm
integer :: flag,k





namelist /parametrosinic/ t_start , t_end
namelist /constantfis/ km, GN, msol, cluz, epsilon0, mu0, hplnck, hbarex, eelect, melect, unosalf, jz0
namelist /parametroscmp/ tol, metodo, orden
!namelist /paramstresh_Newton/ tresh

!open  (unit = 11, file = "datosent/paramstresh_Newton.in", status = "old" )
!read  (unit = 11, *) (tresh(i) , i=1,12)
open  (unit = 10, file = "datosent/paramscmp_Newton.in", status = "old" )
read  (unit = 10, nml = parametroscmp)

open (unit = 22 , file="datosent/parametrosinic-N.in" , status="old")
read (unit = 22 , nml=parametrosinic )
open (unit = 24 , file="datosent/constantesfis.in" , status="old")
read (unit = 24 , nml=constantfis )

close (unit = 10)
close(unit = 11)
close (unit = 22)
close (unit = 24)

open (unit = 101 , file="datossal/salidas-New.f" , status="unknown")

pi = 4.0_wp*atan(1.0_wp)

!write (hache0 ,'( es10 .3) ') h

!hache = 'h = ' // trim ( adjustl ( hache0 ) )

write (* ,*) " "
write (101 ,*) " "

write (* ,*) " metodo =", metodo
write (101 ,*) " metodo =", metodo
write (* ,*) "tol = ",tol
write (101 ,*) "tol = ", tol
write (* ,*) " orden = ", orden
write (101 ,*) " orden = ", orden
write (* ,*) " precicion =", wp
write (101 ,*) " precision =", wp
write (* ,*) " cluz =", cluz
write (101 ,*) "cluz = ", cluz
write (* ,*) " GN =", GN
write (101 ,*) "GN = ", GN
write (* ,*) " msol =", msol
write (101 ,*) "msol = ", msol
write (* ,*) " cluz =", cluz
write (101 ,*) "cluz = ", cluz

unit_len = 1000._wp !km
unit_mass = msol !Masa solar
unit_time = 1._wp !Segundos
Gef = GN * (unit_mass/unit_len**3)  !Constante de la gravedad en km**3/(s**2 * msol)
Gef_c2 = GN * (unit_mass/((cluz**2)*(unit_len**3))) !Constante de la gravedad en km


pi = 4.0_wp*atan(1.0_wp)
M2 = 50._wp
M1 = 40._wp


!Cálculo de los radios de Schwarzschild

r_s1 = Gef_c2*(2*M1)
r_s2 = Gef_c2*(2*M2)


ecc = 0.4_wp
a = 350._wp*r_s1

!precession_deg = (6*pi * Gef_c2 * (M1+M2) / (a*(1-ecc**2))) * (180/pi)


r0 = a * (1.0_wp - ecc)
!v0 = sqrt( (Gef * (M1 + M2) * (1.0_wp + ecc)) / (a * (1.0_wp + ecc)) )

!Cálculo de los períodos orbitales

periodo = 2._wp*pi*sqrt((M1*M2)/Gef*abs(M1*M2)*(M1+M2)) * sqrt(a**3)




y_start(1:6) = (/ (r0* M2)/(M1+M2)  , 0.0_wp, 0.0_wp, &
        -(r0 * M1)/(M1+M2) , 0.0_wp, 0.0_wp /)

y_start(7:12) = &
         (/ 0.0_wp, M2/(M1+M2)*sqrt( Gef*(M1+M2)* (1.0_wp + ecc) / (a * (1.0_wp - ecc) )), 0.0_wp, &
         0.0_wp, -M1/(M1+M2)*sqrt( (Gef*(M1+M2)* (1.0_wp + ecc)) / (a * (1.0_wp - ecc) )) , 0.0_wp /)

dtA0_ini = Delta0_A(y_start(1:3), y_start(4:6), y_start(7:9))
dtA1_ini = Delta1_A(y_start(1:3), y_start(4:6), y_start(7:9), y_start(10:12), M1, M2)
dtAp_ini = DeltaP_A(y_start(1:3), y_start(4:6), y_start(7:9), y_start(10:12), M1, M2)
dtB0_ini = Delta0_B(y_start(4:6), y_start(1:3), y_start(10:12))
dtB1_ini = Delta1_B(y_start(4:6), y_start(1:3), y_start(10:12), y_start(7:9), M1, M2)
dtBp_ini = DeltaP_B(y_start(4:6), y_start(1:3), y_start(10:12), y_start(7:9), M1, M2)

accel_ini_A = Force2_A(y_start(1:3), y_start(7:9), y_start(4:6), y_start(10:12), M1, M2)/M1
accel_ini_B = Force2_B(y_start(4:6), y_start(10:12), y_start(1:3), y_start(7:9), M1, M2)/M2

write (* ,*) " dtA0_ini =" , dtA0_ini
write (101 ,*) "dtA0_ini = " , dtA0_ini

write (* ,*) " dtA1_ini =" , dtA1_ini
write (101 ,*) "dtA1_ini = " , dtA1_ini

write (* ,*) " dtAp_ini =" , dtAp_ini
write (101 ,*) "dtAp_ini = " , dtAp_ini

write (* ,*) " dtB0_ini =" , dtB0_ini
write (101 ,*) "dtB0_ini = " , dtB0_ini

write (* ,*) " dtB1_ini =" , dtB1_ini
write (101 ,*) "dtB1_ini = " , dtB1_ini

write (* ,*) " dtBp_ini =" , dtBp_ini
write (101 ,*) "dtBp_ini = " , dtBp_ini

write (* ,*) " accel_ini_A =" , accel_ini_A
write (101 ,*) "accel_ini_A = " , accel_ini_A

write (* ,*) " accel_ini_B =" , accel_ini_B
write (101 ,*) "accel_ini_B = " , accel_ini_B

write (* ,*) " ecc =" , ecc
write (101 ,*) "ecc = " , ecc

write (* ,*) " a =" , a
write (101 ,*) "a  = " , a

write (* ,*) " r0  =" , r0
write (101 ,*) "r0  = " , r0

!write (* ,*) " v0  =" , v0
!write (101 ,*) "v0  = " , v0

write (* ,*) " M1 (msol) =" , M1
write (101 ,*) "M1 (msol) = " , M1

write (* ,*) " M2 (msol) =" , M2
write (101 ,*) "M2 (msol) = " , M2

write (* ,*) " q1  =" , q1
write (101 ,*) "q1  = " , q1

write (* ,*) " q2  =" , q2
write (101 ,*) "q2  = " , q2

write (* ,*) " Gef_c2 =" , Gef_c2
write (101 ,*) "Gef_c2 = " , Gef_c2

write (* ,*) " periodo =" , periodo
write (101 ,*) "periodo = " , periodo

write (* ,*) " t_end/periodo =" , t_end/(2._wp*pi*sqrt((a**3)*(M1*M2)/((Gef*abs(M1*M2)*(M1+M2)))))
write (101 ,*) "t_end/periodo = " , t_end/(2._wp*pi*sqrt((a**3)*(M1*M2)/((Gef*abs(M1*M2)*(M1+M2)))))

write (* ,*) " r_s1 =" , r_s1
write (101 ,*) "r_s1 = " , r_s1

write (* ,*) " r_s2 =" , r_s2
write (101 ,*) "r_s2 = " , r_s2





write (hache0 ,'( es10 .3) ') h

hache = 'h = ' // trim ( adjustl ( hache0 ) )



write (* ,*) " "
write (101 ,*) " "

write (* ,*) " Gef =" , Gef
write (101 ,*) "Gef = " , Gef

write (* ,*) " Gef_c2 =" , Gef_c2
write (101 ,*) "Gef_c2 = " , Gef_c2

write (* ,*) " r_s1 =" , r_s1
write (101 ,*) "r_s1 = " , r_s1

write (* ,*) " r_s2 =" , r_s2
write (101 ,*) "r_s2 = " , r_s2

write (* ,*) " Precesión teórica:", precession_deg, " grados/orbita"
write (101 ,*) "Precesión teórica:", precession_deg, " grados/orbita"

write (* ,*) " hache0 = " , hache0
write (101 ,*) " hache0 = " , hache0

write (* ,*) " hache = " , hache
write (101 ,*) " hache = " , hache

write (* ,*) " t_start = " , t_start

write (* ,*) " t_end = " , t_end
!write (* ,*) " h = " ,h
write (* ,*) "y_start = " , y_start
write (101 ,*) " t_start = " , t_start
write (101 ,*) " t_end = " , t_end
!write (101 ,*) " h = " ,h
write (101 ,*) " y_start = " , y_start


mcheps = EPSILON (dummy) !It is the smallest value of x such that 1.+x is not equal to 1.
write (*,*)'      mcheps = ',mcheps
dwarf = TINY (dummy)
write (*,*)'       dwarf = ',dwarf !smallest positive number of a numeric data type that can be
                               !represented without loss of precision

write (*,*)   ' sqrt(dwarf) = ',SQRT(dwarf)
write (101,*) '      mcheps = ',mcheps ,"It is the smallest value of x such that 1.+x is not equal to 1."
write (101,*) '       dwarf = ',dwarf , "smallest positive number of a numeric data type that can be"
write (101,*) "                               represented without loss of precision"

write (101,*) ' sqrt(dwarf) = ',SQRT(dwarf)

write (* ,*) " "
write (101 ,*) " "

write (* ,*) " derivt (0._wp , y_start ) = " , derivt (0._wp , y_start )
write (101 ,*) " derivt (0._wp , y_start ) = " , derivt (0._wp , y_start )



write (* ,*) " En0 = " , En0
write (101 ,*) " En0 = " , En0



write (*,*) "Lx_0 = " , Lx_0
write (101,*) "Lx_0 = " , Lx_0
write (*,*) "Ly_0 = " , Ly_0
write (101,*) "Ly_0 = " , Ly_0
write (*,*) "Lz_0 = " , Lz_0
write (101,*) "Lz_0 = " , Lz_0
write (*,*) "Px_0 = " , Px_0
write (101,*) "Px_0 = " , Px_0
write (*,*) "Py_0 = " , Py_0
write (101,*) "Py_0 = " , Py_0
write (*,*) "Pz_0 = " , Pz_0
write (101,*) "Pz_0 = " , Pz_0

close(101)

thres(1:6) = (/ 1e-20_wp, 1e-20_wp, 1e-20_wp, 1e-20_wp, 1e-20_wp, 1e-20_wp /)
thres(7:12) = (/ 1e-20_wp, 1e-20_wp, 1e-20_wp, 1e-20_wp, 1e-20_wp, 1e-20_wp /)
error_assess0 = .false.

t_prev = t_start
archivo = "datossal/Newton_posiciones.dat"
open(unit =10 , file = archivo)

archivo = "datossal/Newton_velocidades.dat"
open(unit =11 , file = archivo)

archivo = 'datossal/New_MA.dat'
open(unit =30 , file = archivo)

archivo = 'datossal/New_momento.dat'
open(unit =40 , file = archivo)

archivo = 'datossal/New_ener.dat'
open(unit =50 , file = archivo)





write(10 ,*) "t", " X_1x " , " X_1y " , " X_1z" , " X_2x " , " X_2y " , " X_2z "
write(11, *) "t", " V_1x ", " V_1y " , " V_1z " , " V_2x " ,"V_2y " , " V_2z "
write (30 ,*) "t" , "Lx", "Ly", "Lz", "(Lz_0 - Lz)/Lz_0"

write (40 ,*) "t" , "Px", "Py", "Pz", "(Py - Py_0)/Py_0"

write (50 ,*) "t" , "En", " (En - En0)/En0 "




call setup(comm,t_start,y_start,t_end,tol,thres,method=metodo,task='step', &
           error_assess=error_assess0)
!           error_assess=.false.)
!           error_assess=.true.)



error_assess0 = .true.

!consistency_done = .false.

do
!   call bsstep(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
    call step_integrate(comm,derivt,t_now,y_now=y_now,flag=flag)

    call constantes_movimiento(y_now, En, Lx, Ly, Lz, Px, Py, Pz, X_CM, Y_CM, Z_CM)
    dy = derivt(t_now,y_now)

    !En_eV = En * (unit_energy_SI / eelect)
    !call open_teval_log()


    write(10,'(E29.8,6(E29.13))') t_now, y_now(1), y_now(2), y_now(3), y_now(4), y_now(5), y_now(6)
    write(11,'(E29.8,6(E29.13))') t_now, y_now(7), y_now(8), y_now(9),y_now(10), y_now(11), y_now(12)
    write(12,'(E29.8,6(E29.13))') t_now, dy(7), dy(8), dy(9),&
                                    dy(10), dy(11), dy(12)
 k = k + 1
   if (t_now==t_end .or. flag > 3) exit
end do
close(10)

close(30)
close(40)
close(50)
close(60)


END PROGRAM newton
