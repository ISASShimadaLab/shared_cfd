#origin=`pwd`
#cd ../../
origin="sample/no13.flame_sheet_chem"

#checkout model
rm checkout.inp
rm checkout_chem.inp
cp $origin/checkout_model.inp .
./checkout.py

#checkout chem
rm checkout_model.inp
cp $origin/checkout_chem.inp .
./checkout.py
mv control_chem.raw.inp          checkout_chem/control_chem.inp
mv chem.inp                      checkout_chem/chem.inp
mv checkout_chem/control.raw.inp checkout_chem/control.inp
cd checkout_chem
make
./driver
