#!/bin/bash

Start_count=0
Stop_count=99
Step_count=1
counter=1
cd ..
mkdir snaps
cd output
for Item in `seq $Start_count $Step_count $Stop_count`; do
     	echo " $Item "
	mkdir test-$Item
	cd test-$Item
	cp ../xsec_rho_$Item rho_input 
	cp ../../src/plot_rxsec_12292020.py .
	python plot_rxsec_12292020.py
	cp pic.png ../../snaps/fig_$Item.png 
	counter=$[$counter +1]
	cd ..
	rm -r test-$Item
done
