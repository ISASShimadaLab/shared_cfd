#origin=`pwd`
#cd ../../
origin="sample/no16.cea_reduction_hp"

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
