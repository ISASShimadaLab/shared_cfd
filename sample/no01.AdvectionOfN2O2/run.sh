#origin=`pwd`
#cd ../../
origin="sample/no01.AdvectionOfN2O2"

# move input file and checkout
cp $origin/*.x .
cp $origin/chem.inp .
cp $origin/checkout.inp .
./checkout.py

# move necessary input files
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
