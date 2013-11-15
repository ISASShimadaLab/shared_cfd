#origin=`pwd`
#cd ../../
origin="sample/no00.shock_tube"

# move input file and checkout
cp $origin/checkout.inp .
./checkout.py

# move necessary input files
cp $origin/control.inp checkout/
rm checkout/control.raw.inp

cp $origin/condition.f90 checkout/
rm checkout/condition.raw.f90

# move to 'checkout' and make/run
cd checkout
make
./main
