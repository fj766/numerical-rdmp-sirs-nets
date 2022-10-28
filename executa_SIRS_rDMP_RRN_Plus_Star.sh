#$ -S /bin/sh

programa=${1}

tam_rede=${2}

RRN_deg=${3}

Star_deg=${4}

lamb0=${5}

divisor=${6}

alp=${7}

ind_lamb=${8}

echo $@ > argumentos.log

time ./${programa} ${tam_rede} ${RRN_deg} ${Star_deg} ${lamb0} ${divisor} ${alp} ${ind_lamb}
