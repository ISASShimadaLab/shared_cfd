directory=debug
#directory=share_code
cp trunk_fs/control.inp       checkout/
cp trunk_fs/control_chem.inp  checkout/
cp trunk_fs/condition.f90     checkout/
cp trunk_fs/restart.bin       checkout/
rm -f  checkout/*raw*
if [ -e checkout/mod_mpi.f90 ]; then
   cp -r checkout ${directory}
   tar czvf share_code.tar.gz ${directory}
   rm -rf ${directory}
   scp share_code.tar.gz maja:/large/y/y535/2dcfd/
fi
