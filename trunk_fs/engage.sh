cp trunk_fs/control.inp       checkout/
cp trunk_fs/control_chem.inp  checkout/
cp trunk_fs/condition.f90     checkout/
cp trunk_fs/restart.bin       checkout/
rm -f  checkout/*raw*
if [ -e checkout/mod_mpi.f90 ]; then
   cp -r checkout share_code
   tar czvf share_code.tar.gz share_code
   rm -rf share_code
   echo "y535 or o384 ('y' or 'o'?)"
   read var
   if [ "${var}" = "o" ]; then
      scp share_code.tar.gz majao:/large/o/o384/
   else
      scp share_code.tar.gz maja:/large/y/y535/2dcfd/
   fi
fi
