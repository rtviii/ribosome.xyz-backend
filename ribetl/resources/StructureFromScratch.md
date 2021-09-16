
# Digitizing the ribosome


A brief note on the repository that will result:


`${ROOT}/static` is the main repository for static files such as structures, individual chains, ligands, tunnel reports etc.

Folder structure is as following:

```typescript
{PDBID}-|
          {PDBID}.json
          {PDBID}.cif
          _ray_{PDBID}.png
          |-CHAINS--|
                    // Individual Subchains
                    |------ {PDBID}_STRAND_{x1}.cif
                    |------ {PDBID}_STRAND_{x2}.cif
                    |------ ...
          |-TUNNEL
                    // Tunnel files
                    |------ {PDBID}_TUNNEL_REPORT.json
                    |------ tunnel{x}.csv
        // Ligand and antibiotics reports
          {PDBID}_LIGAND_{CHEMID1}.json 
          {PDBID}_LIGAND_{CHEMID2}.json 
          ...
```

----------------------

## 1.Semantic Preprocessing

The candidate structures([ex.](candidates_February2020.txt)) are obtained at rcsb from a tabular report on rcsb search api:
https://www.rcsb.org/search, exported as a single txt file with all the structures.
Ideally, there is no preprocessing step to this and we work with the GraphQL endpoint directly and feed it into the database. Far from it currently.

- ```ribxyz``` utility is built to download the relevant fields and transform the rcsb response from a gql, but a bunch of methods are obsoleted. The script is [```generateStructProfileFromRCSBCandidates.sh```](../rcsb-gql-api/requestGqlFrame.ts) is written such that would generate profile given pdbids as lines. *Be careful and execute in root of ribxz*. This produces a template for each structure: the outline of proteins, rnas, ligands. The nomenclature is injencted for proteins where possible.  The resulting templates conform to (our own) schema defined in [ RibosomeTypes.ts](./../../src/../src/RibosomeTypes.ts)

- A bit about the rcsb [ schema ](../src/rcsb-gql-api/archive/expandedRnas-EMDiffnData.gql):

  [ Schema definition ](https://data.rcsb.org/index.html#gql-schema) is provided by RCSB themselves and so far as i can tell, that's what they are using internally too: makes sense to use this here. 
  
```graphql
  entry(entry_id: "3j9m") {
    # Fields relevant to the structure itself.
    rcsb_id
    rcsb_entry_info {
      resolution_combined
    }
    rcsb_external_references {
      link
      type
      id
    }
    exptl {
      method
    }
    # Either diffraction source or 3d-em is present(depending on resolution method)
    diffrn_source {
      details
      pdbx_synchrotron_beamline
      pdbx_synchrotron_site
      pdbx_wavelength
      pdbx_wavelength_list
      source
      type
    }
    em_3d_reconstruction {
      details
      algorithm
      resolution_method
      actual_pixel_size
      id
      resolution
      nominal_pixel_size
      refinement_type
      num_particles
      magnification_calibration
    }
    citation {
      rcsb_authors
      year
      title
      pdbx_database_id_DOI
    }
    struct_keywords {
      pdbx_keywords
      text
    }
    # Proteins and rna : polymer entites
    polymer_entities {
      entry {
        rcsb_id
      }
      pfams {
        rcsb_pfam_accession
        rcsb_pfam_comment
        rcsb_pfam_description
      }
      rcsb_entity_source_organism {
        ncbi_taxonomy_id
        scientific_name
      }
      uniprots {
        rcsb_id
      }
      rcsb_polymer_entity {
        pdbx_description
      }
      # rnas and proteins are differentited here by the  _rcsb_entity_polymer_type_ field
      entity_poly {
        pdbx_seq_one_letter_code
        pdbx_seq_one_letter_code_can
        pdbx_strand_id
        rcsb_entity_polymer_type
        rcsb_sample_sequence_length
        type
      }
    }
    # Ligands and smaller molecules
    nonpolymer_entities {
      pdbx_entity_nonpoly {
        comp_id
        name
        entity_id
      }
      rcsb_nonpolymer_entity {
        formula_weight
        pdbx_description
      }
    }
  }
```

Further, each structure is inducted into the database with the following script: [ _induct_struct.sh ](./../src/resources/inductStructNeo4j.sh). **Check validity of corresponding cypher with the current scheme**. Can of course run in bulk.

### This concludes the semantic part of the preprocessing.


## 2.Structural Preprocessing 

#### These files will be served statically by Django.

Most of these scripts rely on Biopython.PDB's type system and parsers and PyMol's processing capacities.

### 2.1. Generating the structure avatar (MUST)

Ray a pymol image for each struct to serve as an avatar. Would be great to see if a gif can be made: https://pymolwiki.org/index.php/Making_Movies.



### 2.2. Splitting the structure into individual chains (MUST)

The structures aard being [ parsed into individual subchains ](./../ciftools/splitStruct.py ) to be served statically.  
Be wary of structures that pack two models and hence return duplicate chain ids i.e. 1VY5 (AO,CO). These should be split into individual proteins. Not sure yet how to adjust references to them elsewhere in the app.


### 2.3. Chaffing out the ligands where any are present(OPTIONAL)

[ `Ligand.py` ](./../ciftools/Ligand.py) ::parseLigandNeighborhoods


### 2.4. Processing, moving the tunnels into the appropriate folders (OPTIONAL)


This is a curated dataset. Only one of the (typically) 4-20 tunnels produced by MOLE would make the cut and operated on to produce a report of the [ tunnel walls ](./../ciftools/TunnelScripts/WallsReportGeneration.py). This corresponding csv file should also be served statically and needs moving into the correct folder in static.

  The tunnels are created on individual basis. We have some convincing ecoli and thermus thermophilus data and these are the two most represented species in the data, but neither these tunnels perfect nor should we stop here. MOLE did most of the heavy lifting here, but is a tremendous pain to work and interop with, quite hard to configure. 
  Something also has to be done for the structures in which something is present inside the tunnel itself. TUNNELS_LOG.csv 


Ligands must be identified and parsed from those structures which contain them.
___

Given a PDB accession, say __6O8W__ :
- Get the spine of the profile at rcsb: https://www.rcsb.org/pdb/rest/describeMol?structureId=6O8w
RCSB has made some real progress https://data.rcsb.org/index.html#data-schema


# Neo4j/CypherShell Notes

Make sure to __execute neo4j-admin commands as neo4j, the user.__ ```sudo -u neo4j <command>```

:use system

Backup:
neo4j-admin dump --database='...' --to='...date'

*There is a template with all the pfam, uniprot registries and constraints kicking around already*:

  neo4j-admin load --from='_template_rxz.dump' --database=<database> [--force]
