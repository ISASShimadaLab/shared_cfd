#origin=`pwd`
#cd ../../
origin="sample/no12.flame_sheet_flow"

#checkout model
rm checkout.inp
rm checkout_chem.inp
cp $origin/checkout_model.inp .
./checkout.py

#checkout flow
rm checkout_model.inp
cp $origin/checkout.inp .
cp $origin/*.x .
./checkout.py

# move necessary input files
mv control_chem.raw.inp checkout/control_chem.inp
mv chem.inp             checkout/

cp $origin/control.inp checkout/control.inp
rm checkout/control.raw.inp

cp $origin/condition.f90 checkout/
rm checkout/condition.raw.f90

# move to 'checkout' and make/run
cd checkout
make
ulimit -s unlimited
./main
