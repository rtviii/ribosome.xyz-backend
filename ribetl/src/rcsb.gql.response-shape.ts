export interface Polymer_Entity {
  entry:{
    rcsb_id:string
  }

  rcsb_polymer_entity_container_identifiers: {
        asym_ids     : string[]
        auth_asym_ids: string[]
        entry_id     : string
        entity_id    : string
      }

  pfams: {
    rcsb_pfam_accession  : string;
    rcsb_pfam_comment    : string;
    rcsb_pfam_description: string;
    }[] | null;

  rcsb_entity_source_organism: { 
    ncbi_taxonomy_id: number;
    scientific_name : string;
    }[];

  uniprots: {rcsb_id: string;}[] | null;

  rcsb_polymer_entity: {
    pdbx_description: string;
  };

  entity_poly: {
    pdbx_seq_one_letter_code    : string;
    pdbx_seq_one_letter_code_can: string;
    pdbx_strand_id              : string;
    rcsb_entity_polymer_type    : string;
    rcsb_sample_sequence_length : number;
    type                        : string;
  };
}
export interface Nonpolymer_Entity {

  pdbx_entity_nonpoly: {
    comp_id  : string;
    name     : string;
    entity_id: string;
  };



  rcsb_nonpolymer_entity: {
    formula_weight          : number;
    pdbx_description        : string;
    pdbx_number_of_molecules: number;
    details                 : string
  };

}

export interface PDBGQLResponse {
  entry: {
    rcsb_id: string;
    struct_keywords: {
      pdbx_keywords: string
      text         : string
    };
    rcsb_entry_info: { resolution_combined: number[]; };

    rcsb_external_references: {link:string; type:string; id:string}[]

    exptl                   : { method: string }[];

    citation: {
      rcsb_authors        : string[];
      year                : number;
      title               : string;
      pdbx_database_id_DOI: string;
    }[];

    polymer_entities   : Polymer_Entity[];
    nonpolymer_entities: Nonpolymer_Entity[] | null;
  };
}
