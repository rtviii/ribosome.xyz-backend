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
        struct.rcsb_external_ref_id                    = CASE WHEN ref_id                = null then \"null\" else ref_id END,
        struct.rcsb_external_ref_type                  = CASE WHEN ref_type              = null then \"null\" else ref_type END,
        struct.rcsb_external_ref_link                  = CASE WHEN ref_link              = null then \"null\" else ref_link END" | cypher-shell 