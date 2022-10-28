#!/bin/bash

modulos='mod_types.f90 mod_rndgen_multiple.f90 geraRede.f90 mod_tools_redes.f90'

principal='sirs_message_passing_estacionario.f90'

exe='mp_10k_g28_1'

#flags='-traceback -check all -heap-arrays -O3 -g -fp-stack-check -check uninit -check bounds -ftrapuv'
 
flags=''

ifort ${flags} ${modulos} ${principal} -o ${exe}

chmod +x ${exe}

rm -r *.mod &> erroRM.log

time ./${exe}
