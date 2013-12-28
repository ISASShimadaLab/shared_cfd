#origin=`pwd`
#cd ../../
origin="sample/no23.chemkin_mono_uv"
#checkout chem
rm checkout.inp
rm checkout_model.inp
cp $origin/chem.inp .
cp $origin/checkout_chem.inp .
./checkout.py
cp $origin/control_chem.inp checkout_chem/
rm checkout_chem/control_chem.raw.inp
cp $origin/control.chem.inp checkout_chem/control.inp
rm checkout_chem/control.raw.inp
cd checkout_chem
make
./driver

