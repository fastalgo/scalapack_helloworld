rm execute
source /opt/intel/bin/compilervars.sh intel64
#source /opt/intel/impi/5.0.2.044/bin64/mpivars.sh
source /opt/intel/impi_latest/bin64/mpivars.sh
CC -o execute pdgesv.cpp -mkl:cluster
sbatch test.sl
