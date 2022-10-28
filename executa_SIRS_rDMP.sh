#$ -S /bin/sh

programa=${1}
tam_rede=${2}
grau_RRN=${3}
lamb0=${4}
divisor=${5}
lambdaf=${6}
alp=${7}
ind_lamb=${8}
ind_soCalculaIPR=${9}


echo $@ > argumentos.log

time ./${programa} ${tam_rede} ${grau_RRN} ${lamb0} ${divisor} ${lambdaf} ${alp} ${ind_lamb} ${ind_soCalculaIPR}
