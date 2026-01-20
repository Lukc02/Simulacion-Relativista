
rm -f *.o *.mod modulos/*.o modulos/*.mod Tesis/*.o Tesis/*.mod

gfortran -g -fcheck=all -Ofast -fimplicit-none\
         modulos/modulo-precision.f90\
         modulos/modulo-constantes-newton.f90\
         modulos/modulo-variables.f90\
         modulos/Funciones.f90\
         modulos/fuerza_retardada_lorentz.f90\
         modulos/fuerza_retardada_orden2.f90\
         modulos/modulo-lagran4.f90\
         modulos/modulo-retardo.f90\
         modulos/modulo_derivada_lorentz.f90\
         modulos/rksuite_90_r1.f90\
         Tesis/Newton_rk_suite.f90\
         -o Newton.x

echo
