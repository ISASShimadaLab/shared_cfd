#origin=`pwd`
#cd ../../
origin="sample/no02.boundary_layer"

# move input file and checkout
cp $origin/checkout.inp .
cp $origin/*.x .
./checkout.py

# move necessary input files
cp $origin/chem.inp checkout/

cp $origin/control_chem.inp checkout/
rm checkout/control_chem.raw.inp

cp $origin/control.inp checkout/
rm checkout/control.raw.inp

cp $origin/condition.f90 checkout/
rm checkout/condition.raw.f90

# move to 'checkout' and make/run
cd checkout
make
ulimit -s unlimited
./main
