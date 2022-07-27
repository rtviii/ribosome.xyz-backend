#!/usr/bin/bash
set -e
# set -u
set -o pipefail

# For remote access:
# echo "match (n) return count(n)" | cypher-shell -a neo4j+s://ribosome.xyz:7687 --format plain -u 'rt' -p 'rrr' --database 'riboauth'


usage (){
      echo "

      _____________

      The following env vars must be set:
      NEO4J_CURENTDB=<password>
      NEO4J_URI=<server>

      ______________
      
      $0  /
      [-f] [path to a structure profile(.json)] /
      [-a] [Neo4j remote server i.e. localhost:7678] /
      [-d] [Database name i.e. 'ribolocal'] /
      " >&2
      exit 1
}

structure() {
  if [ -f $FILEPATH ];
  then
    file=$(basename -- $FILEPATH)
    extension=${file##*.}
    if [ $extension != "json" ];
    then
      echo "The profile file must be a .json. Exiting."
      exit $((2))
    fi
    structid=${file::4}
    structid=${structid^^}
  else
    echo "$FILEPATH is not an acceptable file"
    exit $((1))
  fi


  echo "Commiting structure $structid core to $NEO4J_URI."
  echo "call apoc.load.json(\"file:///static/$structid/$file\") yield value
  with                                                       value.rcsb_id as pdbid,
                                                            value.expMethod as exp,
                                                            value.resolution as reso,

                                                            value.rcsb_external_ref_id   as ref_id  ,
                                                            value.rcsb_external_ref_type as ref_type,
                                                            value.rcsb_external_ref_link as ref_link,

                                                            value.citation_year as cit_year,
                                                            value.citation_rcsb_authors as cit_authors,
                                                            value.citation_title as cit_title,
                                                            value.citation_pdbx_doi as cit_doi,

                                                            value.pdbx_keywords_text as kwordstext,
                                                            value.pdbx_keywords as kwords,

                                                            value.src_organism_ids   as src_organism_ids  ,
                                                            value.src_organism_names as src_organism_names,

                                                            value.host_organism_ids     as host_organism_ids    ,
                                                            value.host_organism_names   as host_organism_names  ,


                                                            value

  merge                 ( struct                                 :RibosomeStructure{
          rcsb_id               : pdbid                                  ,
          expMethod             : exp                                    ,
          resolution            : reso                                   ,
          citation_rcsb_authors : cit_authors                            ,
          citation_title        : cit_title                              ,

          pdbx_keywords     : kwords     ,
          pdbx_keywords_text: kwordstext,

          src_organism_ids           : src_organism_ids                                  ,
          src_organism_names         : src_organism_names                                ,

          host_organism_ids           : host_organism_ids                                  ,
          host_organism_names         : host_organism_names

          })

          on                      create                                  set
          struct.rcsb_external_ref_id                    = CASE WHEN ref_id                = null then \"null\" else ref_id END,
          struct.rcsb_external_ref_type                  = CASE WHEN ref_type              = null then \"null\" else ref_type END,
          struct.rcsb_external_ref_link                  = CASE WHEN ref_link              = null then \"null\" else ref_link END,
          struct.citation_pdbx_doi                       = CASE WHEN cit_doi              = null then \"null\" else ref_link END,
          struct.citation_year                          = CASE WHEN cit_year              = null then \"null\" else cit_year END" | \
          cypher-shell -a "$NEO4J_URI" --format plain -u $NEO4J_USER -p $NEO4J_PASSWORD --database $NEO4J_CURENTDB

}

rnas (){
    
  if [ -f $FILEPATH ];
  then
    file=$(basename $FILEPATH)
    extension=${file: -4}
    if [ $extension != "json" ];
    then
      echo "The profile file must be a .json. Exiting."
      exit 2
    fi
    structid=${file::4}
    structid=${structid^^}
  else
    echo "$FILEPATH is not an acceptable file"
    exit -1
  fi


  echo "Commiting RNA of $structid to $NEO4J_URI."
  echo "call apoc.load.json(\"file:///static/$structid/$file\") yield value
  with value 
      unwind                                 value .rnas as rna
      merge                               (  newrna:RNA {

      asym_ids                         : rna.asym_ids,
      auth_asym_id                    : rna.auth_asym_id,

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
  merge (n)-[:belongs_to]-(nc)" | \
  cypher-shell -a "$NEO4J_URI" --format plain -u $NEO4J_USER -p $NEO4J_PASSWORD --database $NEO4J_CURENTDB




}

proteins(){

  if [ -f $FILEPATH ];
  then
    file=$(basename -- $FILEPATH)
    extension=${file##*.}
    if [ $extension != "json" ];
    then
      echo "The profile file must be a .json. Exiting."
      exit $((2))
    fi
    structid=${file::4}
    structid=${structid^^}
  else
    echo "$FILEPATH is not an acceptable file"
    exit $((1))
  fi


  echo "Commiting proteins of $structid to $NEO4J_URI."
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
  merge (n)-[:member_of]->(nc)" | \
  cypher-shell -a "$NEO4J_URI" --format plain -u $NEO4J_USER -p $NEO4J_PASSWORD --database $NEO4J_CURENTDB
}

ligands (){
  if [ -f $FILEPATH ];
  then
    file=$(basename -- $FILEPATH)
    extension=${file##*.}
    if [ $extension != "json" ];
    then
      echo "The profile file must be a .json. Exiting."
      exit $((2))
    fi
    structid=${file::4}
    structid=${structid^^}
  else
    echo "$FILEPATH is not an acceptable file"
    exit $((1))
  fi

  echo "Commiting ligands of $structid to $NEO4J_URI."
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
  return s.rcsb_id, newligand.chemicalId " | \
  cypher-shell -a "$NEO4J_URI" --format plain -u $NEO4J_USER -p $NEO4J_PASSWORD --database $NEO4J_CURENTDB

}


check_if_exists (){
  ifexists_cypher="match (n:RibosomeStructure{rcsb_id:\"$1\"}) return (n)"
  response=$( echo $ifexists_cypher | cypher-shell -a "$NEO4J_URI" --format plain -u $NEO4J_USER -p $NEO4J_PASSWORD --database $NEO4J_CURENTDB )
  if [[ ${#response} -gt 1 ]] ; then
    echo "Something with the rcsb_id $1 already exists in the database. Failing on purpose.."
    echo "Used query: $ifexists_cypher"
    echo $ifexists_cypher
    exit 1
  fi


}


DEFAULT_NEO4J_ADDRESS="127.0.0.1:7687"
while getopts 'd:f:a:' OPTION; do
  case "$OPTION" in
    d)
      NEO4J_CURENTDB="$OPTARG"
      ;;
    f)
      FILEPATH="$OPTARG"
      if [ -f $FILEPATH ];
      then
        file=$FILEPATH
        extension=${file##*.}
        if [ $extension != "json" ];
        then
          echo "The profile file must be a .json. Exiting."
          exit $((2))
        fi
      else
        echo "$FILEPATH is not an acceptable file"
        exit $((1))
      fi
      echo "Got filepath $FILEPATH"
      ;;
    a)
      NEO4J_URI="$OPTARG"
      ;;
    ?)
    usage
      ;;
  esac
done


if [[ -z $FILEPATH ]] || [[ -z $NEO4J_CURENTDB ]] || [[ -z $NEO4J_URI ]] 
then 
  usage 
else 

  if [ -f $FILEPATH ];
  then
    file=$(basename -- $FILEPATH)
    extension=${file##*.}
    if [ $extension != "json" ];
    then
      echo "The profile file must be a .json. Exiting."
      exit $((2))
    fi
    structid=${file::4}
    structid=${structid^^}
  else
    echo "$FILEPATH is not an acceptable file"
    exit $((1))
  fi

  check_if_exists $structid;
  structure;
  proteins;
  rnas;
  ligands;
fi
