


# Database Template

Refer to [ `neo4j-steps.md` ](./neo4j-steps.mdj) for database creation.

The data that constitute the backbone of the database are found in `/ribetl/resources/cumulativeData`. This should be present in the neo4j's import folder when creating the database. Typically `/var/lib/neo4j/import`

Create the base-types, constraints of the database with scripts contained in [ this cypher ](../ribetl/resources/cypher-tools/current.modular.cypher).


# ETL


	When performing operations inside the neo4j directories, ex. import, act as neo4j the user to avoid permissions issues: sudo su neo4j

process the download in bulk by applying driver.ts in parallel to all the pdbidi in the download file: eliminate the commas with ```sed -i "s/\,/\\n/g" rcsb_pdb_ids_20210926175604.txt```

## Graph Profiles 

ENS are in  [ src/driver.ts ](../ribetl/src/driver.ts)

```typescript
 dotenv.config({ path: '/home/rxz/dev/ribetl/.env' });
```


The API's response is transformed by the [ driver scripts ](../ribetl/src/requestGQLProfile.ts) to the form appropriate to the database injestion. For a given `$RCSBID` The result is stored in `/static/$RCSBID/$RCSBID.json`.4


***Propagate any changes to the types to the INDUCT SCRIPTS, to the FRONT-END interface***.
1. We really value well-defined interfaces and the correspondence of types. Hence, they are specified in [ ```RibosomeTypes.ts``` ](../ribetl/src/RibosomeTypes.ts). This is the basis for the Neo4j ontology as well as a guiding structure for the front-end's datatypes. **If changes/additions are to be made changes/additions to the application, the ought to begin here.** 

2. RCSB GraphQL endpoint is at `https://data.rcsb.org/graphql`. It is queried with the desired shaped for each molecule. 
Template query is in [ template_query ](../ribetl/src/requestGQLProfile.ts). The response shape should conform to the **types**(see 1.). 

3. The resultant `.json` profile is used to initiate nodes and links in the database. For a given structure, the script creates individual components in sequence. Refer to [ cypher ](../ribetl/resources/cypher-tools/current.modular.cypher).

	Currently and roughly:

		1. (Merge)Create the structure node if one doesn't exist
		2. Find this struct, for each  contained protein -- create its node and connect to struct
		3. Connect proteins to PFAM Families
		4. Connect proteins to nomenclature classes
		5. Find this struct, for each  contained rna -- create its node and connect to struct
		6. Connect rnas to nomenclature classes
		7. Connect ligands


There is some ambiguity right now as to what to consider a Ligand. Some elongation factors land in the RP category because of their classification as a polypeptide.


## Structural Files

Actual RCSB structures are stored in [ `batch_download` ](../ribetl/batch_download/) folder (update on each cycle).  Query is of the form: [all ribosome structures, of resolution smaller than 4 angstrom deposited 2014.](https://www.rcsb.org/search?request=%7B%22query%22%3A%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22struct_keywords.pdbx_keywords%22%2C%22operator%22%3A%22contains_phrase%22%2C%22negation%22%3Afalse%2C%22value%22%3A%22RIBOSOME%22%7D%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%2C%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_entry_info.resolution_combined%22%2C%22operator%22%3A%22less_or_equal%22%2C%22negation%22%3Afalse%2C%22value%22%3A4%7D%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%2C%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_accession_info.initial_release_date%22%2C%22operator%22%3A%22greater%22%2C%22negation%22%3Afalse%2C%22value%22%3A%222014-01-01T00%3A00%3A00Z%22%7D%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%5D%2C%22label%22%3A%22text%22%7D%5D%7D%2C%22return_type%22%3A%22entry%22%2C%22request_info%22%3A%7B%22query_id%22%3A%220a8b586b4227c759d60304b8272bb0d3%22%7D%2C%22request_options%22%3A%7B%22pager%22%3A%7B%22start%22%3A0%2C%22rows%22%3A25%7D%2C%22scoring_strategy%22%3A%22combined%22%2C%22sort%22%3A%5B%7B%22sort_by%22%3A%22score%22%2C%22direction%22%3A%22desc%22%7D%5D%7D%7D)

1. Each one has to have its chains renamed according to the new Ban nomenclature (assigned during graph profile generation). [ Scripts here ](../ribetl/ciftools/renaming_structs/). To be deposited at `/static/$RCSBID/$RCSBID.cif`

2. Process the binding sites of ligands, elongation factors etc.  (including whatever else gets included). [ Scripts ](../ribetl/ciftools/binding_site.py). To be deposited at `/static/$RCSBID/LIGAND_$LIGANDID*.json`

3. Split the structure into individual protien and rna to be deposited at `/static/$RCSBID/CHAINS/$RCSBID_STRAND_$STRANID*.cif`.


# Todos

- [ ] Mitochondrial category
- [ ] Stages of the translation cycle



```zsh
find . -name "*.json" | xargs grep 'rcsb_pdbx_description'  | awk -F  ':' ' $3 !~ /protein |RNA|rRNA|PROTEIN|Protein|mS|uL|UL|eL|bL|bS|BS|uS|eS|bS|mL|EL|rna|protein|ul|ml|RACK1/  {print $3}'
```