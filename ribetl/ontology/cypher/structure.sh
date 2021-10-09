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
        struct.citation_year                          = CASE WHEN cit_year              = null then \"null\" else cit_year END" | cypher-shell  --database riboxyz
