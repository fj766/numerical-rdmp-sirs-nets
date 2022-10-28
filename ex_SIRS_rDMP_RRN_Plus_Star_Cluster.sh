#!/bin/bash

ind_lamb=${1}

num_base=1

tam_rede=${num_base}0000

RRN_deg=6

Star_deg=998

lamb0=0.0

divisor=5.0

lambdaf=1.0

alp=0.5

echo $PWD
#=======================================================================
#maquina='cluster'
maquina='casa'
#=======================================================================
#flags='-check all -traceback'

#flags='-heap-arrays -O3 -g -fp-stack-check'

#flags='-traceback -check all -heap-arrays -O3 -g -fp-stack-check'

flags=''
#=======================================================================
nzeros=$( echo "scale=0; l(${tam_rede}/${num_base})/l(10) " | bc -l)

IFS='.'

arr=( ${gama_exp} )

gama_sp=${arr[0]}${arr[1]}

alp_arr=( ${alp} )

IFS=' '
##########################################################
dependencias='mod_rndgen_multiple.f90 geraRede.f90 mod_tools_redes.f90 mod_types.f90'
principal='main_SIRS_rDMP_Simplificado_RRN_Plus_Star.f90'

executavel='rDMP_RRN_Plus_Star_n'${num_base}${nzeros}'RRN_deg_'${RRN_deg}'Star_deg'${Star_deg}'a'${alp_arr[0]}'p'${alp_arr[1]}'_indLambd_'${ind_lamb}

exec_cluster='executa_SIRS_rDMP_RRN_Plus_Star.sh'

rm ${executavel} > erroRmExe.log

ifort ${dependencias} ${principal} ${flags} -o ${executavel}

rm -r *.mod

#=======================================================================
if [[ ${maquina} == 'casa' ]];then
   #--------------------------------------------------------------------
   time ./${executavel} ${tam_rede} ${RRN_deg} ${Star_deg} ${lamb0} ${divisor} ${alp} ${ind_lamb}
   #--------------------------------------------------------------------
elif [[ ${maquina} == 'cluster' ]];then
   #--------------------------------------------------------------------
   qsub -N ${executavel} -cwd ${exec_cluster} ${executavel} ${tam_rede} ${RRN_deg} ${Star_deg} ${lamb0} ${divisor} ${alp} ${ind_lamb}
   #--------------------------------------------------------------------
else
   #--------------------------------------------------------------------
   echo "Maquina nao-reconhecida. Abortando script..."
   #--------------------------------------------------------------------
   exit
fi
