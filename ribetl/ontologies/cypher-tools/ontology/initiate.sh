#!/usr/bin/bash

# Make sure to move 'resources' repository from ribetl to neo4j import folder
# (that's where file:// reads from by default : /var/lib/neo4j/import)

# Create constraints
echo '
CREATE CONSTRAINT IF NOT EXISTS ON (ipro:InterProFamily   ) ASSERT ipro.family_id  IS UNIQUE;
CREATE CONSTRAINT IF NOT EXISTS ON (go  :GOClass          ) ASSERT go.class_id   IS UNIQUE;
CREATE CONSTRAINT IF NOT EXISTS ON (pf  :PFAMFamily       ) assert pf.family_id  is unique;

CREATE CONSTRAINT IF NOT EXISTS ON (q   :RibosomeStructure) Assert q.rcsb_id    IS UNIQUE;

CREATE CONSTRAINT IF NOT EXISTS ON (pc  :ProteinClass ) assert pc.class_id   is unique;
CREATE CONSTRAINT IF NOT EXISTS ON (rc  :RNAClass     ) assert rc.class_id   is unique;
CREATE CONSTRAINT IF NOT EXISTS ON (lig :Ligand       ) assert lig.chemicalId is unique;

' | cypher-shell;


# INTERPRO families
echo 'call apoc.load.json("file:///resources/cumulativeData/interpro.json") yield value
with value as v
merge (q:InterProFamily{ family_id:KEYS(v)[0],type:v[KEYS(v)[0]].type,description:v[KEYS(v)[0]].name})' | cypher-shell


# INTERPRO-GeneOntology connections
for i in 1 2 3 4;
do 
	echo "call apoc.load.json(\"file:///resources/cumulativeData/interpro-go/part$i.json\") yield value as go
	merge (inode:InterProFamily{family_id:go.InterPro})
	merge (gonode:GOClass{go_class:go.GO})
	on create set gonode.annotation = go.GO_annotation
	merge (inode)-[:interpro_go{annotation:go.interpro_class}]-(gonode)" | cypher-shell;
done


# PFAM-Interpro connections
for r in 1 2 3 4;
do 
	echo "call apoc.load.json(\"file:///resources/cumulativeData/pfam-interpro/part$r.json\") yield value as entry
	with entry.metadata as datum
	with datum where datum.integrated is not null
	merge (inode:InterProFamily{family_id: datum.integrated})
	merge (pnode:PFAMFamily{family_id: datum.accession, family_type:datum.type})
	merge (inode)-[:interpro_pfam]-(pnode);" | cypher-shell;
done


# Small Subunit nomenclature map
echo "call apoc.load.json(\"file:///resources/SSUMap.json\") yield value
unwind(keys(value)) as key
merge (pc:ProteinClass {class_id:key});" | cypher-shell;

# Large Subunit nomenclature map
echo "\
call apoc.load.json(\"file:///resources/LSUMap.json\") yield value
unwind(keys(value)) as key
merge (pc:ProteinClass {class_id:key})" | cypher-shell;

# RNA nomenclature classes
echo 'UNWIND [ 
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
  merge (n:RNAClass {class_id:rnaclass})' | cypher-shell;
