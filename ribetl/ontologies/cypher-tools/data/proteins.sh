#!/usr/bin/bash

NEOIMPORT='/var/lib/neo4j/import'

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
with value                              , struct
     unwind                              value.proteins as protein
with protein                            ,
     value                              ,
     struct
     merge                               (rp:RibosomalProtein {
     asym_ids:      protein.asym_ids,
     auth_asym_ids:protein.auth_asym_ids,
        
     parent_rcsb_id                      : protein.parent_rcsb_id,

     pfam_comments                       : protein.pfam_comments,
     pfam_descriptions                   : protein.pfam_descriptions,
     pfam_accessions                     : protein.pfam_accessions,

     rcsb_source_organism_description    : protein.rcsb_source_organism_description,
     rcsb_source_organism_id             : protein.rcsb_source_organism_id,
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
with rp, value
match(s:RibosomeStructure {rcsb_id: value.rcsb_id})
create (protein)-[:RibosomalProtein_of]->(s)

with rp,struct,value
     unwind  rp    .   pfam_accessions as pfamils
     match  (pf    :   PFAMFamily      {family_id:pfamils})
with rp,struct,value,pf
     merge  (rp    )-[:Belogns_To     ]->(pf)


match (n:RibosomalProtein) where n.nomenclature[0] is not null
merge (nc:RPClass{class_id:n.nomenclature[0]})
merge (n)-[:BelongsTo]-(nc)" | cypher-shell
