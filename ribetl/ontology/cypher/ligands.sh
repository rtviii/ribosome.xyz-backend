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


echo "Looking at structure $structid"
echo "call apoc.load.json(\"file:///static/$structid/$file\") yield value
with value.rcsb_id as struct, value
       unwind           value.ligands as lig
       merge            (newligand:Ligand {
	chemicalId          : lig.chemicalId         ,
	chemicalName        : lig.chemicalName       ,
	formula_weight      : lig.formula_weight     ,
	pdbx_description    : lig.pdbx_description   ,
	number_of_instances : lig.number_of_instances
       })

  with newligand, value
  match (s:RibosomeStructure {rcsb_id: value.rcsb_id})
  merge            (newligand)<-[:contains_ligand]-(s)
  return s.rcsb_id, newligand.chemicalId " | cypher-shell --database riboxyz
