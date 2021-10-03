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


# For               every structure, point the script to a json profile. Induct its elements sequentially:
# RibosomeStructure node
# Protein           nodes
# RNA               nodes
# Ligand            nodes

#Throwing struct into cypher-shell
echo "call apoc.load.json(\"file:///static/$structid/$file\") yield value
with                                                       value.rcsb_id as pdbid,
                                                           value.expMethod as exp,
                                                           value.resolution as reso,

                                                           value.rcsb_external_ref_id as ref_id,
                                                           value.rcsb_external_ref_type as ref_type,
                                                           value.rcsb_external_ref_link as ref_link,

                                                           value.citation_year as cit_year,
                                                           value.citation_rcsb_authors as cit_authors,
                                                           value.citation_title as cit_title,
                                                           value.citation_pdbx_doi as cit_doi,

                                                           value.pdbx_keywords_text as kwordstext,
                                                           value.pdbx_keywords as kwords, 

                                                           value._organismId as orgid,
                                                           value._organismName as orgname,
                                                           value

merge                 ( struct                                 :RibosomeStructure{
        rcsb_id               : pdbid                                  ,
        expMethod             : exp                                    ,
        resolution            : reso                                   ,
        citation_year         : cit_year                               ,
        citation_rcsb_authors : cit_authors                            ,
        citation_title        : cit_title                              ,
        citation_pdbx_doi     : cit_doi                                ,

        pdbx_keywords     : kwords     ,
        pdbx_keywords_text: kwordstext,

        _organismId           : orgid                                  ,
        _organismName         : orgname                                
        
        })

        on                      create                                  set
        struct .  rcsb_external_ref_id                    = CASE WHEN ref_id                = null then "null" else ref_id END,
        struct .  rcsb_external_ref_type                  = CASE WHEN ref_type              = null then "null" else ref_type END,
        struct .  rcsb_external_ref_link                  = CASE WHEN ref_link              = null then "null" else ref_link END
with value                              , struct
     unwind                              value.proteins as protein
with protein                            ,
     value                              ,
     struct
     merge                               (rp:RibosomalProtein {
        //      Check proteins
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

})-[:RibosomalProtein_of                ]->(struct)
    on create                set
    rp.rcsb_pdbx_description = CASE WHEN protein.rcsb_pdbx_description = null then "null" else protein.rcsb_pdbx_description END
//     connect protein to PFAMFamily
with rp    , struct,   value
     unwind  rp    .   pfam_accessions as pfamils
     match  (pf    :   PFAMFamily      {family_id:pfamils})
with rp    , struct,   value          ,pf
     merge  (rp    )-[:Belogns_To     ]->(pf)

// Connect RNAS
with value                              ,   struct
     unwind                                 value .rnas as rna
     merge                               (  newrna:rRNA {
     asym_ids                         : rna.asym_ids,
     auth_asym_ids                    : rna.auth_asym_ids,
     parent_rcsb_id                      :  rna   .parent_rcsb_id,

     rcsb_source_organism_description    :  rna   .rcsb_source_organism_description,
     rcsb_source_organism_id             :  rna   .rcsb_source_organism_id,
     nomenclature                        :  rna.nomenclature,

     entity_poly_strand_id               :  rna   .entity_poly_strand_id,
     entity_poly_seq_one_letter_code     :  rna   .entity_poly_seq_one_letter_code,
     entity_poly_seq_one_letter_code_can:   rna   .entity_poly_seq_one_letter_code_can,
     entity_poly_seq_length              :  rna   .entity_poly_seq_length,
     entity_poly_polymer_type            :  rna   .entity_poly_polymer_type,
     entity_poly_entity_type             :  rna   .entity_poly_entity_type
})-[:rRNA_of                            ]->(struct)
     on                                     create set newrna.rcsb_pdbx_description = CASE WHEN rna.rcsb_pdbx_description = null then "null" else rna.rcsb_pdbx_description END
// connect Ligands
with   value           , struct
       unwind           value.ligands as lig
       merge            (l:Ligand {
       chemicalId      : lig.chemicalId,
       chemicalName    : lig.chemicalName,
       formula_weight  : lig.formula_weight,
       pdbx_description: lig.pdbx_description})
       merge            (l)<-[:ContainsLigand]-(struct)
return {struct: struct.rcsb_id}
" | cypher-shell -u rt -p rrr  --format plain
