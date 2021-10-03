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
with value 
     unwind                                 value .rnas as rna
     merge                               (  newrna:RNA {

     asym_ids                         : rna.asym_ids,
     auth_asym_ids                    : rna.auth_asym_ids,

     parent_rcsb_id:  rna.parent_rcsb_id,

     nomenclature:  rna.nomenclature,
     
     ligand_like: rna.ligand_like,

     src_organism_ids  :rna.src_organism_ids,
     src_organism_names:rna.src_organism_names,
     host_organism_ids    :rna.host_organism_ids  ,
     host_organism_names  :rna.host_organism_names,


     entity_poly_strand_id               :  rna.entity_poly_strand_id,
     entity_poly_seq_one_letter_code     :  rna.entity_poly_seq_one_letter_code,
     entity_poly_seq_one_letter_code_can :  rna.entity_poly_seq_one_letter_code_can,
     entity_poly_seq_length              :  rna.entity_poly_seq_length,
     entity_poly_polymer_type            :  rna.entity_poly_polymer_type,
     entity_poly_entity_type             :  rna.entity_poly_entity_type
}) on create set newrna.rcsb_pdbx_description = CASE WHEN rna.rcsb_pdbx_description = null then \"null\" else rna.rcsb_pdbx_description END

with newrna, value
match(s:RibosomeStructure {rcsb_id: value.rcsb_id})
create (newrna)-[:rna_of]->(s);

match (n:RNA) where n.nomenclature[0] is not null
merge (nc:RNAClass{class_id:n.nomenclature[0]})
merge (n)-[:belongs_to]-(nc)" | cypher-shell 