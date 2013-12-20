#origin=`pwd`
#cd ../../
origin="sample/no22.cea_flow"

#checkout model
rm checkout.inp
rm checkout_chem.inp
cp $origin/checkout_model.inp .
./checkout.py

#checkout chem
rm checkout_model.inp
cp $origin/checkout_chem.inp .
./checkout.py
cp $origin/control_chem.inp checkout_chem/
rm control_chem.raw.inp
cp $origin/control.reduction.inp checkout_chem/control.inp
rm checkout_chem/control.raw.inp
cd checkout_chem
make
./driver
cp chem.inp.new ../chem.inp
cd ..

#checkout flow
rm checkout_chem.inp
cp $origin/checkout.inp .
cp $origin/*.x .
./checkout.py

# move necessary input files
cp $origin/control_chem.inp checkout/

cp $origin/control.flow.inp checkout/control.inp
rm checkout/control.raw.inp

cp $origin/condition.f90 checkout/
rm checkout/condition.raw.f90

# move to 'checkout' and make/run
cd checkout
make
ulimit -s unlimited
./main
