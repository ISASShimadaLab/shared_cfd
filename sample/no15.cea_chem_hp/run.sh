#origin=`pwd`
#cd ../../
origin="sample/no15.cea_chem_hp"

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
mv checkout_chem/control.raw.inp checkout_chem/control.inp
cd checkout_chem
make
./driver
