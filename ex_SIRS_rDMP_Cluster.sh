#!/bin/bash

num_base=1

tam_rede=${num_base}0000

grau_RRN=6

lamb0=0.0

divisor=5.0

lambdaf=1.0

alp=0.5
ind_lamb=${1}

ind_soCalculaIPR=0
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
principal='main_SIRS_rDMP_Simplificado.f90'

executavel='rDMP_n'${num_base}${nzeros}g${gama_sp}'a'${alp_arr[0]}'p'${alp_arr[1]}'_ams_'${ind_amostra}'_indLambd_'${ind_lamb}


rm ${executavel} > erroRmExe.log

ifort ${dependencias} ${principal} ${flags} -o ${executavel}

rm -r *.mod

#=======================================================================
if [[ ${maquina} == 'casa' ]];then
   #--------------------------------------------------------------------
   time ./${executavel} ${tam_rede} ${grau_RRN} ${lamb0} ${divisor} ${lambdaf} ${alp} ${ind_lamb} ${ind_soCalculaIPR}
   #--------------------------------------------------------------------
elif [[ ${maquina} == 'cluster' ]];then
   #--------------------------------------------------------------------
   qsub -N ${executavel} -cwd executa_SIRS_rDMP.sh ${executavel} ${tam_rede} ${grau_RRN} ${lamb0} ${divisor} ${lambdaf} ${alp} ${ind_lamb} ${ind_soCalculaIPR}
   #--------------------------------------------------------------------
else
   #--------------------------------------------------------------------
   echo "Maquina nao-reconhecida. Abortando script..."
   #--------------------------------------------------------------------
   exit
fi
