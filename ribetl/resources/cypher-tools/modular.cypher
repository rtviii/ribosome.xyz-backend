CREATE CONSTRAINT ON (ipro:InterProFamily) ASSERT ipro.family_id IS UNIQUE
CREATE CONSTRAINT ON (go:GOClass) ASSERT go.class_id IS UNIQUE
CREATE CONSTRAINT ON (q:RibosomeStructure) Assert q.rcsb_id IS UNIQUE
CREATE CONSTRAINT ON (pf:PFAMFamily) assert pf.family_id is unique
CREATE CONSTRAINT ON (lig:Ligand) assert lig.chemicalId is unique
CREATE CONSTRAINT ON (nc:NomenclatureClass) assert nc.class_id is unique

call apoc.load.json("file:///resources/cumulativeData/interpro.json") yield value
with value as v
merge (q:InterProFamily{ family_id:KEYS(v)[0],type:v[KEYS(v)[0]].type,description:v[KEYS(v)[0]].name})

CALL apoc.load.json('file:///resources/cumulativeData/interpro-go/part{1-4}.json') yield value as go
merge (inode:InterProFamily{family_id:go.InterPro})
merge (gonode:GOClass{go_class:go.GO})
on create set gonode.annotation = go.GO_annotation
merge (inode)-[:mp_InterPro_GO{annotation:go.interpro_class}]-(gonode)

call apoc.load.json("file:///resources/cumulativeData/pfam-interpro/part{1-4}.json") yield value as entry
with entry.metadata as datum
with datum where datum.integrated is not null
merge (inode:InterProFamily{family_id: datum.integrated})
merge (pnode:PFAMFamily{family_id: datum.accession, family_type:datum.type})
merge (inode)-[:mp_InterPro_PFAM]-(pnode)

call apoc.load.json('file:///resources/SSUMap.json') yield value
unwind(keys(value)) as key
merge (nc:NomenclatureClass {class_id:key})

call apoc.load.json('file:///resources/LSUMap.json') yield value
unwind(keys(value)) as key
merge (nc:NomenclatureClass {class_id:key})

UNWIND [ 
  "5SrRNA"  ,
  "5.8SrRNA",
  "12SrRNA" ,
  "16SrRNA" ,
  "21SrRNA" ,
  "23SrRNA" ,
  "25SrRNA" ,
  "28SrRNA" ,
  "35SrRNA" ,
  "mRNA"    ,
  "tRNA"    ]  as rnaclass
  create (n:NomenclatureClass {class_id:rnaclass})

// following the type
call apoc.load.json("file:///static/$structid/$file") yield value
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
with struct


call apoc.load.json("file:///static/$structid/$file") yield value
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
})
on create                set
rp.rcsb_pdbx_description = CASE WHEN protein.rcsb_pdbx_description = null then "null" else protein.rcsb_pdbx_description END
with rp, value
match(s:RibosomeStructure {rcsb_id: value.rcsb_id})
create (protein)-[:RibosomalProtein_of]->(s)

//     connect protein to PFAMFamily
with rp,struct,value
     unwind  rp    .   pfam_accessions as pfamils
     match  (pf    :   PFAMFamily      {family_id:pfamils})
with rp,struct,value,pf
     merge  (rp    )-[:Belogns_To     ]->(pf)
// Connect RNAS


call apoc.load.json("file:///static/7OF4/7OF4.json") yield value
with value 
     unwind                                 value .rnas as rna
     Merge                               (  newrna:rRNA {
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
})on create set newrna.rcsb_pdbx_description = CASE WHEN rna.rcsb_pdbx_description = null then "null" else rna.rcsb_pdbx_description END
with newrna, value
match(s:RibosomeStructure {rcsb_id: value.rcsb_id})
create (newrna)-[:rRNA_of                            ]->(s)



// connect Ligands
call apoc.load.json("file:///static/7OF4/7OF4.json") yield value
       unwind           value.ligands as lig
       merge            (l:Ligand {
       chemicalId      : lig.chemicalId,
       chemicalName    : lig.chemicalName,
       formula_weight  : lig.formula_weight,
       pdbx_description: lig.pdbx_description})

with l, value
match(s:RibosomeStructure {rcsb_id: value.rcsb_id})
create            (l)<-[:ContainsLigand]-(s)
return struct;

// ---------------------------------------------

// CONNECT NOMENCLATURE
match (n:RibosomalProtein) where n.nomenclature[0] is not null
merge (nc:RPClass{class_id:n.nomenclature[0]})
merge (n)-[:BelongsTo]-(nc)

// CONNECT NOMENCLATURE
match (n:rRNA) where n.nomenclature[0] is not null
merge (nc:RNAClass{class_id:n.nomenclature[0]})
merge (n)-[:BelongsTo]-(nc)






