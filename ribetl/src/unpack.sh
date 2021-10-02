#!/usr/bin/zsh


parallel  'chmod +x {1}; unzip {1}' ::: $(ls *.zip)
parallel  'gunzip {1}; echo "Processed {}"' ::: $(ls *.gz)


for i in $( ls *.cif) 
do
	RCSBID=${i:0:4:u};
	mkdir $RCSBID; mv $i $RCSBID/$RCSBID.cif;
done
