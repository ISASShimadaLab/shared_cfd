name=debug
all:
clean:
	rm -rf checkout checkout*.inp *.x *.dat run.sh chem.inp control*.inp *.bak.* grid_bak checkout_chem *.pyc
pack:
	cp -r checkout $(name)
	tar czvf share_code.tar.gz $(name)
	rm -rf $(name)
	scp share_code.tar.gz maja:/large/y/y535/2dcfd/
