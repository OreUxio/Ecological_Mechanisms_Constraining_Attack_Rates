#!/bin/sh
> list_of_files_run.txt
mkdir Data/Raw/Full_assembly_model/Old
mv Data/Raw/Full_assembly_model/* Data/Raw/Full_assembly_model/
for ((a=1;a<10;a++)); do
	for b in 4; do
        	for c in 250; do
			target=$a
			target+=_
			target+=$b
			target+=_
			target+=$c
			targets=$target
			targets+=.cfg
			mkdir .../Data/$target
			cp NewAsMini.cfg important/Data1/$target/$targets
			sed -i "/totally_random_link_strength_sigma/s/= .*/= $b/" important/Data1/$target/$targets
			sed -i "/target_S/s/= .*/= $c/" important/Data1/$target/$targets
			echo $target >> list_of_files_run.txt
		done
	done
done
