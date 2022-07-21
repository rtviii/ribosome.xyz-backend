#!/usr/bin/bash

if [ $# -lt $((1))  ];
then 
        echo "Not enough arguments!"
        exit $((1))
fi

filepath=$1
if [ -f $filepath ];
then
	file=$(basename -- $filepath)
	extension=${file##*.}
	if [ $extension != "json" ];
	then
		echo "The profile file must be a .json. Exiting."
		exit $((2))
	fi
	structid=${file::4}
	structid=${structid^^}
     echo "Processing $structid"
else
	echo "$filepath is not an acceptable file"
	exit $((1))
fi


echo "call apoc.load.json(\"file:///static/$structid/$file\") yield value

with value, value.rcsb_id as struct
     unwind                              value.proteins as protein
with protein                            ,
     value                              ,
     struct
     merge                               (rp:Protein {
     asym_ids:      protein.asym_ids,
     auth_asym_id:protein.auth_asym_id,
        
     parent_rcsb_id                      : protein.parent_rcsb_id,

     pfam_comments                       : protein.pfam_comments,
     pfam_descriptions                   : protein.pfam_descriptions,
     pfam_accessions                     : protein.pfam_accessions,

     src_organism_ids  :protein.src_organism_ids,
     src_organism_names:protein.src_organism_names,
     host_organism_ids    :protein.host_organism_ids  ,
     host_organism_names  :protein.host_organism_names,
     
     ligand_like:protein.ligand_like,

     uniprot_accession                   : protein.uniprot_accession,

     rcsb_pdbx_description              : protein.rcsb_pdbx_description,

     entity_poly_strand_id              : protein.entity_poly_strand_id,
     entity_poly_seq_one_letter_code    : protein.entity_poly_seq_one_letter_code,
     entity_poly_seq_one_letter_code_can: protein.entity_poly_seq_one_letter_code_can,
     entity_poly_seq_length             : protein.entity_poly_seq_length,
     entity_poly_polymer_type           : protein.entity_poly_polymer_type,
     entity_poly_entity_type            : protein.entity_poly_entity_type,

     nomenclature                        : protein.nomenclature
})
on create                set
rp.rcsb_pdbx_description = CASE WHEN protein.rcsb_pdbx_description = null then \"null\" else protein.rcsb_pdbx_description END

with rp, value, struct

match(s:RibosomeStructure {rcsb_id: value.rcsb_id})
create (rp)-[:protein_of]->(s)
with rp,struct,value
     unwind  rp    .   pfam_accessions as pfamils
     match  (pf    :   PFAMFamily      {family_id:pfamils})

with rp,struct,value,pf
     merge  (rp    )-[:belongs_to     ]->(pf);

match (n:Protein) where n.nomenclature[0] is not null
merge (nc:ProteinClass {class_id:n.nomenclature[0]})
merge (n)-[:member_of]->(nc)" | cypher-shell --database $RIBOXYZ_DB_NAME
# You must set the database name in the script : "export RIBOXYZ_DB_NAME= XXX"

