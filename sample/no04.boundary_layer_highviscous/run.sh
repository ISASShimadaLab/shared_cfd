#origin=`pwd`
#cd ../../
origin="sample/no04.boundary_layer_highviscous"

# move input file and checkout
cp $origin/checkout.inp .
cp $origin/*.x .
./checkout.py

# move necessary input files
cp $origin/control.inp checkout/
rm checkout/control.raw.inp

cp $origin/condition.f90 checkout/
rm checkout/condition.raw.f90

cp $origin/thermal_model.f90 checkout/

# move to 'checkout' and make/run
cd checkout
make
ulimit -s unlimited
./main
