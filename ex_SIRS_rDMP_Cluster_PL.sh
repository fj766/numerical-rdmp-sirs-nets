#!/bin/bash

ind_amostra=${1}
ind_lamb=${2}

ind_noh=$((${ind_lamb} % 10))
num_base=1

tam_rede=${num_base}000000

grau_min=3; gama_exp=3.5

lamb0=0.0

divisor=5.0

lambdaf=1.0

alp=0.5

ind_soCalculaIPR=0
#=======================================================================
#maquina='cluster'
#=======================================================================
maquina='cluster_ind'
maq=4
noh=(0-0 0-1 0-12 0-14 0-2 0-20 0-3 0-5 0-6 0-8)
no=${noh[${ind_noh}]}
echo ${no}
#=======================================================================
#maquina='casa'
#=======================================================================
#flags='-check all -traceback'

#flags='-heap-arrays -O3 -g -fp-stack-check'

#flags='-traceback -check all -heap-arrays -O3 -ipo -g -fp-stack-check'

flags='-Ofast -ipo'
#=======================================================================
nzeros=$( echo "scale=0; l(${tam_rede}/${num_base})/l(10) " | bc -l)

IFS='.'

arr=( ${gama_exp} )

gama_sp=${arr[0]}${arr[1]}

alp_arr=( ${alp} )

IFS=' '
##########################################################
dependencias='mod_rndgen_multiple.f90 geraRede.f90 mod_tools_redes.f90 mod_types.f90'
principal='main_SIRS_rDMP_Simplificado_Y4_Teorico.f90'

executavel='rDMP_n'${num_base}${nzeros}g${gama_sp}'a'${alp_arr[0]}'p'${alp_arr[1]}'_ams_'${ind_amostra}'_indLambd_'${ind_lamb}


rm ${executavel} > erroRmExe.log

ifort ${dependencias} ${principal} ${flags} -o ${executavel}

rm -r *.mod

#=======================================================================
if [[ ${maquina} == 'casa' ]];then
   #--------------------------------------------------------------------
   time ./${executavel} ${ind_amostra} ${tam_rede} ${grau_min} ${gama_exp} ${lamb0} ${divisor} ${lambdaf} ${alp} ${ind_lamb} ${ind_soCalculaIPR}
   #--------------------------------------------------------------------
elif [[ ${maquina} == 'cluster' ]];then
   #--------------------------------------------------------------------
   qsub -N ${executavel} -cwd executa_SIRS_rDMP_PL.sh ${executavel} ${ind_amostra} ${tam_rede} ${grau_min} ${gama_exp} ${lamb0} ${divisor} ${lambdaf} ${alp} ${ind_lamb} ${ind_soCalculaIPR}
   #--------------------------------------------------------------------
elif [[ ${maquina} == 'cluster_ind' ]];then
   qsub -N ${executavel} -q all.q@compute-${no} -cwd executa_SIRS_rDMP_PL.sh ${executavel} ${ind_amostra} ${tam_rede} ${grau_min} ${gama_exp} ${lamb0} ${divisor} ${lambdaf} ${alp} ${ind_lamb} ${ind_soCalculaIPR}
   echo "Escolhemos o noh "${no}
else
   #--------------------------------------------------------------------
   echo "Maquina nao-reconhecida. Abortando script..."
   #--------------------------------------------------------------------
   exit
fi
