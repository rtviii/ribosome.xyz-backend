import axios from "axios";
import { uniq } from "lodash";
import {
  BanClass,
  Ligand,
  RibosomalProtein,
  RibosomeStructure,
  RNAClass,
  rRNA,
} from "./RibosomeTypesv2";
import { small_subunit_map } from "./../resources/small-subunit-map";
import { large_subunit_map } from './../resources/large-subunit-map';
import {
  Nonpolymer_Entity,
  PDBGQLResponse,
  Polymer_Entity,
} from "./rcsb.gql.response-shape";



const query_template = (pdbid: string) => {
  pdbid = pdbid.toUpperCase();
  return `https://data.rcsb.org/graphql?query=
  {
  entry(entry_id: "${pdbid}") {
    rcsb_id
    struct_keywords {
      pdbx_keywords
      text
    }
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
    polymer_entities {
      entry{
        rcsb_id
      }
      rcsb_polymer_entity_container_identifiers {
        asym_ids
        auth_asym_ids
        entry_id
        entity_id
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
      entity_poly {
        pdbx_seq_one_letter_code
        pdbx_seq_one_letter_code_can
        pdbx_strand_id
        rcsb_entity_polymer_type
        rcsb_sample_sequence_length
        type
      }
    }
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
      rcsb_nonpolymer_entity_container_identifiers {
        chem_ref_def_id
        nonpolymer_comp_id
        prd_id
        rcsb_id
      }
    }
  }
}`;
};

const matchRPNomenclature = (
  polymer: Polymer_Entity
): BanClass[] => {

  var banregex = /([ueb][ls]\d{1,2})/gi
  // * check authors's annotations. if classes are present --> use that.
  var finds = banregex.exec(polymer.rcsb_polymer_entity.pdbx_description)
  if (finds !== null) {
    var firstcap = finds[0]
    var name = (firstcap[0].toLowerCase() + firstcap[1].toUpperCase() + firstcap.slice(2)) as BanClass
    return [name]
  }
  // * then resort to pfams
  if (!polymer.pfams) {
    return []
  } else {
    var pfamIds = polymer.pfams.reduce<string[]>((acc, item) => {
      acc.push(item.rcsb_pfam_accession);
      return acc;
    }, []);
    var nomenclature: BanClass[] = [];
    pfamIds!.map(id => {
      Object.entries(large_subunit_map).map(entry => {
        if (entry[1].pfamDomainAccession.includes(id)) { nomenclature.push(entry[0] as BanClass); }
      });
      Object.entries(small_subunit_map).map(entry => {
        if (entry[1].pfamDomainAccession.includes(id)) {
          nomenclature.push(entry[0] as BanClass);
        }
      });
    });

    return uniq(nomenclature)
  }
};



const matchRNANomenclature = (polymer: Polymer_Entity): RNAClass[] => {


  const rna_reg: Record<RNAClass, RegExp> = {
    "5SrRNA": /\b(5s)/gi,
    "5.8SrRNA": /\b(5\.8s)/gi,
    "12SrRNA" : /\b(12s)/gi,
    "16SrRNA" : /\b(16s)/gi,
    "21SrRNA" : /\b(21s)/gi,
    "23SrRNA" : /\b(23s)/gi,
    "25SrRNA" : /\b(25s)/gi,
    "28SrRNA" : /\b(28s)/gi,
    "35SrRNA" : /\b(35s)/gi,
    "mRNA"    : /(mrna)|\b(messenger)\b/gi,
    "tRNA"    : /(trna)|\b(transfer)\b/gi,
    }
  
  var rnatypes = Object.keys(rna_reg) as RNAClass[]
  for (var c of rnatypes){
    var matches = rna_reg[c].exec(polymer.rcsb_polymer_entity.pdbx_description)
    if (matches){
      return [ c as RNAClass ] 
    }
  }
  return []

}


const inferOrganismFromPolymers: (
  proteins: RibosomalProtein[]
) => [string[], number[]] = proteins => {
  var structOrgIds  : number[] = [];
  var structOrgNames: string[] = [];

  proteins.map(protein => {
    protein.rcsb_source_organism_description
      ? structOrgNames.push(...protein.rcsb_source_organism_description)
      : "";
    protein.rcsb_source_organism_id
      ? structOrgIds.push(...protein.rcsb_source_organism_id)
      : "";
  });
  return [uniq(structOrgNames), uniq(structOrgIds)];
};

const extractRefs = (
  external_refs: { link: string; type: string; id: string }[]
) => {
  var externalRefIds  : string[] = [];
  var externalRefTypes: string[] = [];
  var externalRefLinks: string[] = [];

  external_refs.map(ref => {
    externalRefIds.push(ref.id);
    externalRefTypes.push(ref.type);
    externalRefLinks.push(ref.link);
  });
  return [externalRefIds, externalRefTypes, externalRefLinks];
};

const reshape_ToLigand = (nonpoly: Nonpolymer_Entity): Ligand => {
  return {
    pdbx_description    : nonpoly.rcsb_nonpolymer_entity.pdbx_description        ,
    formula_weight      : nonpoly.rcsb_nonpolymer_entity.formula_weight          ,
    chemicalId          : nonpoly.pdbx_entity_nonpoly   .comp_id                 ,
    chemicalName        : nonpoly.pdbx_entity_nonpoly   .name                    ,
    number_of_instances : nonpoly.rcsb_nonpolymer_entity.pdbx_number_of_molecules
  };
};

const reshape_PolyEntity_to_rRNA =(plm:Polymer_Entity):rRNA =>{

      var organism_ids          : number[] ;
      var organism_descriptions : string[] ;

  if (plm.rcsb_entity_source_organism) {

      organism_ids = plm.rcsb_entity_source_organism.map(org => org.ncbi_taxonomy_id);
      organism_descriptions = plm.rcsb_entity_source_organism.map(org => org.scientific_name)
  }else{

 organism_ids = []; organism_descriptions = [] 
  }
    

  return {
    nomenclature                       : matchRNANomenclature(plm),
    asym_ids                           : plm                  .rcsb_polymer_entity_container_identifiers.asym_ids,
    auth_asym_ids                      : plm                  .rcsb_polymer_entity_container_identifiers.auth_asym_ids,
    parent_rcsb_id                     : plm                  .entry.rcsb_id,
    rcsb_source_organism_id            : organism_ids,
    rcsb_source_organism_description   : organism_descriptions,
    rcsb_pdbx_description              : plm                  .rcsb_polymer_entity.pdbx_description,
    entity_poly_strand_id              : plm                  .entity_poly.pdbx_strand_id,
    entity_poly_seq_one_letter_code    : plm                  .entity_poly.pdbx_seq_one_letter_code,
    entity_poly_seq_one_letter_code_can: plm                  .entity_poly.pdbx_seq_one_letter_code_can,
    entity_poly_seq_length             : plm                  .entity_poly.rcsb_sample_sequence_length,
    entity_poly_entity_type            : plm                  .entity_poly.type,
    entity_poly_polymer_type           : plm                  .entity_poly.rcsb_entity_polymer_type
  }
}

const reshape_PolyEntity_to_RibosomalProtein = (plm:Polymer_Entity):RibosomalProtein =>{

      if (plm.pfams) {

        var pfam_comments = uniq(
          plm.pfams.map(pfam => pfam.rcsb_pfam_comment)
        );
        var pfam_descriptions = uniq(
          plm.pfams.map(pfam => pfam.rcsb_pfam_description)
        );
        var pfam_accessions = uniq(
          plm.pfams.map(pfam => pfam.rcsb_pfam_accession)
        );

      } 
      else {
        var pfam_comments    : string[] = [];
        var pfam_descriptions: string[] = [];
        var pfam_accessions  : string[] = [];
      }

      var organism_ids         : number[];
      var organism_descriptions: string[];

      organism_ids          = plm.rcsb_entity_source_organism.map(org => org.ncbi_taxonomy_id);
      organism_descriptions = plm.rcsb_entity_source_organism.map(org => org.scientific_name);

  return {

    nomenclature                       : matchRPNomenclature(plm),
    asym_ids                           : plm                  .rcsb_polymer_entity_container_identifiers.asym_ids,
    auth_asym_ids                      : plm                  .rcsb_polymer_entity_container_identifiers.auth_asym_ids,
    parent_rcsb_id                     : plm                  .entry.rcsb_id,
    pfam_accessions                    : pfam_accessions,
    pfam_comments                      : pfam_comments,
    pfam_descriptions                  : pfam_descriptions,
    rcsb_source_organism_id            : organism_ids,
    rcsb_source_organism_description   : organism_descriptions,
    uniprot_accession                  : plm                  .uniprots ? plm.uniprots.map(entry => entry.rcsb_id)    : [],
    rcsb_pdbx_description              : plm                  .rcsb_polymer_entity.pdbx_description,
    entity_poly_strand_id              : plm                  .entity_poly.pdbx_strand_id,
    entity_poly_seq_one_letter_code    : plm                  .entity_poly.pdbx_seq_one_letter_code,
    entity_poly_seq_one_letter_code_can: plm                  .entity_poly.pdbx_seq_one_letter_code_can,
    entity_poly_seq_length             : plm                  .entity_poly.rcsb_sample_sequence_length,
    entity_poly_entity_type            : plm                  .entity_poly.type,
    entity_poly_polymer_type           : plm                  .entity_poly.rcsb_entity_polymer_type
  }
}

export const processPDBRecord = async (
  pdbid: string
): Promise<RibosomeStructure> => {

  return await axios.get(query_template(pdbid))
  .then(response => {
    var pdbRecord: PDBGQLResponse = response .data .data                                                                                       ;
    var            proteins       = pdbRecord.entry.polymer_entities   .filter(poly => poly.entity_poly.rcsb_entity_polymer_type === "Protein");
    var            rnas           = pdbRecord.entry.polymer_entities   .filter(poly => poly.entity_poly.rcsb_entity_polymer_type === "RNA")    ;
    var            ligands        = pdbRecord.entry.nonpolymer_entities                                                                        ;

    var reshaped_proteins : RibosomalProtein []        = proteins                         .map(protein => reshape_PolyEntity_to_RibosomalProtein (protein ));
    var reshaped_rrnas    : rRNA             []        = rnas                             .map(rna     => reshape_PolyEntity_to_rRNA             (rna     ));
    var reshaped_ligands  : Ligand           [] | null = ligands  == null ? null : ligands.map(r       => reshape_ToLigand                       (r       ));

    var [organismNames, organismIds] = inferOrganismFromPolymers(reshaped_proteins);
    var externalRefs = extractRefs(pdbRecord.entry.rcsb_external_references);
    var pub          = pdbRecord.entry.citation[0];

    var kwords_text  = pdbRecord.entry.struct_keywords
      ? pdbRecord.entry.struct_keywords.text
      : null;
    var kwords  = pdbRecord.entry.struct_keywords
      ? pdbRecord.entry.struct_keywords.pdbx_keywords
      : null;

    var reshaped: RibosomeStructure = {
      rcsb_id               : pdbRecord.entry.rcsb_id,
      expMethod             : pdbRecord.entry.exptl[0].method,
      resolution            : pdbRecord.entry.rcsb_entry_info.resolution_combined[0],

      rcsb_external_ref_id  : externalRefs[0],
      rcsb_external_ref_type: externalRefs[1],
      rcsb_external_ref_link: externalRefs[2],

      citation_year        : pub.year,
      citation_rcsb_authors: pub.rcsb_authors,
      citation_title       : pub.title,
      citation_pdbx_doi    : pub.pdbx_database_id_DOI,

      pdbx_keywords_text: kwords_text,
      pdbx_keywords     : kwords,
      _organismId       : organismIds,
      _organismName     : organismNames,

      proteins: reshaped_proteins,
      rnas    : reshaped_rrnas,
      ligands : reshaped_ligands,
    };

    return reshaped;
  });

};
