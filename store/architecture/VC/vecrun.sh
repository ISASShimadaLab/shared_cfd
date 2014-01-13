#!/bin/sh
qsub <<!QSUB
#PBS -q D
#PBS -l cpunum_job=8
#PBS -l memsz_job=3.0gb
#PBS -l elapstim_req=600
#PBS -mbe
#!/bin/sh

cd /large/y/y535/share_code
./main
!QSUB
#PBS -v F_FTRACE=YES
