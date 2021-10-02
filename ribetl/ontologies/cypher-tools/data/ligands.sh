#!/usr/bin/bash


filepath=$1
if [ -f $filepath ];
then
	file=$(basename $filepath)
	extension=${file: -4}
	if [ $extension != "json" ];
	then
		echo "The profile file must be a .json. Exiting."
		exit 2
	fi
	structid=${file::4}
	structid=${structid^^}
else
	echo "$filepath is not an acceptable file"
	exit -1
fi


echo "call apoc.load.json(\"file:///static/$structid/$file\") yield value
with value.rcsb_id as struct, value

       unwind           value.ligands as lig
       merge            (l:Ligand {
       chemicalId          : lig.chemicalId         ,
       chemicalName        : lig.chemicalName       ,
       formula_weight      : lig.formula_weight     ,
       pdbx_description    : lig.pdbx_description   ,
       number_of_instances : lig.number_of_instances
       })
with l, value,struct
match(s:RibosomeStructure {rcsb_id: value.rcsb_id})
create            (l)<-[:contains_ligand]-(s)
return l.chemicalId, s.rcsb_id;" | cypher-shell 
