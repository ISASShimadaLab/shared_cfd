#!/bin/sh
. /opt/JXlocal/bin/rjob_set

qsub <<!QSUB
#@\$-q QJOB
#@\$-r share_code
#@\$-lP Nproc
#@\$-lm 5gb
#@\$-cp 600
#@\$-x
#@\$-mb -me
#@\$
#!/bin/sh

cd \${QSUB_WORKDIR}

FLIB_FASTOMP=TRUE ;export FLIB_FASTOMP

mpiexec -n Nproc ./main
RC=\$?

if [ \$RC -gt 0 ]; then
        echo "return code: \$RC"
fi

exit \$RC

!QSUB

#fpcoll -Icall,balance,hwm,mpi -l30 -o result.txt mpiexec -n Nproc ./main
