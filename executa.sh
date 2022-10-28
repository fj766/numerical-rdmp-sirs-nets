#$ -S /bin/sh

programa=${1}

echo $@ > argumentos.log

time ./${programa}
