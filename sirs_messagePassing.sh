#!/bin/bash

modulos='mod_types.f90 mod_rndgen_multiple.f90 geraRede.f90 mod_numerico.f90'

principal='sirs_message_passing.f90'

exe='sirs_message_passing'

flags= #'-traceback -check all -heap-arrays -O3 -g -fp-stack-check'


ifort ${flags} ${modulos} ${principal} -o ${exe} #&> erroCompilacao.log

chmod +x ${exe}

rm -r *.mod &> erroRM.log
