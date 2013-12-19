#origin=`pwd`
#cd ../../
origin="sample/no17.cea_plot_hp"

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
cp $origin/control.inp checkout_chem/
rm checkout_chem/control.raw.inp
cd checkout_chem
make
./driver
