#origin=`pwd`
#cd ../../
origin="sample/no12.flame_sheet_chem"

#rm
rm checkout.inp
cp $origin/checkout_chem.inp .
./checkout.py

mv checkout_chem/control_chem.raw.inp checkout_chem/control_chem.inp
mv checkout_chem/control.raw.inp      checkout_chem/control.inp
cd checkout_chem
make
./driver
cd ..

# move input file and checkout
cp checkout_chem/chem.inp .
cp $origin/checkout.inp .
cp $origin/*.x .
./checkout.py

# move necessary input files
cp checkout_chem/control_chem.inp checkout/

cp $origin/control.inp checkout/control.inp
rm checkout/control.raw.inp

cp $origin/condition.f90 checkout/
rm checkout/condition.raw.f90

# move to 'checkout' and make/run
cd checkout
make
ulimit -s unlimited
./main
