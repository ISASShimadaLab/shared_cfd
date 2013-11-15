origin=`pwd`
cd ../../

# move input file and checkout
cp $origin/checkout.inp .
cp $origin/debug0.x .
./checkout.py

# move necessary input files
cp $origin/chem.inp checkout/

cp $origin/control.inp checkout/
rm checkout/control.raw.inp

cp $origin/control_chem.inp checkout/
rm checkout/control_chem.raw.inp

cp $origin/condition.f90 checkout/
rm checkout/condition.raw.f90

# move to 'checkout' and make/run
cd checkout
make
./main
