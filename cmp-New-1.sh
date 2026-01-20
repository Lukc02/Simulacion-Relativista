
#rm -f *.o *.mod modulos/*.o modulos/*.mod Tesis/*.o Tesis/*.mod

gfortran -g -fcheck=all -Ofast -fimplicit-none\
         modulos/modulo-precision.f90\
         modulos/modulo-constantes-newton.f90\
         modulos/modulo-variables.f90\
         modulos/Funciones.f90\
         modulos/fuerza_retardada_lorentz.f90\
         modulos/fuerza_gravitacional.f90\
         modulos/modulo-retardo-newton.f90\
         modulos/Derivada.f90\
         modulos/rksuite_90_r1.f90\
         Tesis/Newton.f90\
         -o Newton-1.x

echo
