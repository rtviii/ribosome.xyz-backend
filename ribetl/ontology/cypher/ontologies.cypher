call apoc.load.json("file:///resources/cumulativeData/interpro.json") yield value
with value as v
merge (q:InterProFamily{ family_id:KEYS(v)[0],type:v[KEYS(v)[0]].type,description:v[KEYS(v)[0]].name})


CALL apoc.load.json('file:///resources/cumulativeData/interpro-go/part{1-4}.json') yield value as go
merge (inode:InterProFamily{family_id:go.InterPro})
merge (gonode:GOClass{go_class:go.GO})
on create set gonode.annotation = go.GO_annotation
merge (inode)-[:mp_InterPro_GO{annotation:go.interpro_class}]-(gonode);

call apoc.load.json("file:///resources/cumulativeData/pfam-interpro/part{1-4}.json") yield value as entry
with entry.metadata as datum
with datum where datum.integrated is not null
merge (inode:InterProFamily{family_id: datum.integrated})
merge (pnode:PFAMFamily{family_id: datum.accession, family_type:datum.type})
merge (inode)-[:mp_InterPro_PFAM]-(pnode);

call apoc.load.json('file:///resources/SSUMap.json') yield value
unwind(keys(value)) as key
merge (nc:NomenclatureClass {class_id:key});

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
