from operator import sub
from numpy import log
from rest_framework.decorators import api_view
from rest_framework.response import Response
from neo4j import  Result, GraphDatabase
import os
#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯
uri        =  os.getenv( 'NEO4J_URI'                                      )
authglobal = (os.getenv( 'NEO4J_USER'      ),os.getenv( 'NEO4J_PASSWORD' ))
current_db =  os.getenv( 'NEO4J_CURRENTDB'                                )

#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯

def _neoget(CYPHER_STRING:str):
    driver = GraphDatabase.driver(uri, auth= authglobal )
    def parametrized_query(tx, **kwargs):
        result:Result = tx.run(CYPHER_STRING, **kwargs)
        return result.value()

    with driver.session() as session:
        return session.read_transaction(parametrized_query)

#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯

#? ---------------------------LIGANDS
@api_view(['GET', 'POST'])
def get_individual_ligand(request):
    params = dict(request.GET)
    chemId = params['chemId'][0]
    cypher ="""
    match (n:Ligand{{chemicalId:"{}"}}) return n;""".format(chemId)
    return Response(_neoget(cypher))


@api_view(['GET', 'POST'])
def get_all_ligands(request):

    CYPHER_STRING="""
        match (l:Ligand)-[]-(r:RibosomeStructure)  where 
        not l.chemicalName  contains "ION" 
        and not l.chemicalName contains "CLUSTER"
        and not l.chemicalName contains "["
        return {{  
        polymer    : false,
        description: l.chemicalName,
        chemicalId : l.chemicalId,
        presentIn  : {{
                src_organism_ids: r.src_organism_ids,
                description     : l.chemicalName,
                citation_title  : r.citation_title,
                expMethod       : r.expMethod,
                rcsb_id         : r.rcsb_id,
                resolution      : r.resolution
            }}
        }}
    """.format()
    return Response(_neoget(CYPHER_STRING))


@api_view(['GET', 'POST'])
def get_all_ligandlike(request):
    CYPHER_STRING = """
        match (l {{ligand_like:true}})-[]-(r:RibosomeStructure) 
        return {{
            polymer     : true,
            description : l.rcsb_pdbx_description,
            presentIn  : {{
                auth_asym_id    : l.auth_asym_id,
                src_organism_ids: r.src_organism_ids,
                description     : l.rcsb_pdbx_description,
                citation_title  : r.citation_title,
                expMethod       : r.expMethod,
                rcsb_id         : r.rcsb_id,
                resolution      : r.resolution
            }}
        }}""".format()
    return Response(_neoget(CYPHER_STRING))

#? ---------------------------STRUCTS



# -+=-=-=-=-=-=-=-=-

@api_view(['GET'])
def get_RibosomeStructure(request): #<------------- This ought to return an object that conforms to a Ribosome Structure type defined in redux/types
    params = dict(request.GET)
    pdbid  = str.upper(params['pdbid'][0])
    cypher = """
    match (n:RibosomeStructure{{rcsb_id:"{pdbid}"}})
    optional match (rr:RNA)-[]-(n)
    with n, collect(rr) as rrna
    optional match (rp:Protein)-[]-(n)
    with n, rrna,  collect(rp) as rps
    optional match (l:Ligand)-[]-(n)
    with n, rrna, rps, collect(l) as ligs
    return {{
        
                expMethod             : n.expMethod             ,
                resolution            : n.resolution            ,

                pdbx_keywords         : n.pdbx_keywords         ,
                pdbx_keywords_text    : n.pdbx_keywords_text    ,

                rcsb_external_ref_id  : n.rcsb_external_ref_id  ,
                rcsb_external_ref_type: n.rcsb_external_ref_type,
                rcsb_external_ref_link: n.rcsb_external_ref_link,

                citation_year         : n.citation_year         ,
                citation_rcsb_authors : n.citation_rcsb_authors ,
                citation_title        : n.citation_title        ,
                citation_pdbx_doi     : n.citation_pdbx_doi     ,

                src_organism_ids      : n.src_organism_ids      ,
                src_organism_names    : n.src_organism_names    ,

                host_organism_ids     : n.host_organism_ids     ,
                host_organism_names   : n.host_organism_names   ,

                proteins : rps,
                rnas     : rrna,
                ligands  : ligs}}
        
        

    """.format_map({"pdbid":pdbid})
    return Response(_neoget(cypher))

@api_view(['GET']) 
def get_ligands_by_struct(request):
    CYPHER_STRING="""match (n:RibosomeStructure)-[]-(l:Ligand)
           return {{ title: n.citation_title, struct:n.rcsb_id, organism:n.src_organism_names, taxid:n.src_organism_ids, 
           ligands:collect({{ chemid: l.chemicalId, name:l.chemicalName, number:number_of_instances }})}}""".format_map({})
    return Response(_neoget(CYPHER_STRING))

@api_view(['GET'])
def get_struct(request):
    params = dict(request.GET)
    pdbid  = str.upper(params['pdbid'][0])
    cypher = """
    match (n:RibosomeStructure{{rcsb_id:"{pdbid}"}})
    optional match (rr:RNA)-[]-(n)
    with n, collect(rr) as rrna
    optional match (rp:Protein)-[]-(n)
    with n, rrna,  collect(rp) as rps
    optional match (l:Ligand)-[]-(n)
    with n, rrna, rps, collect(l) as ligs
    return {{structure: n, ligands: ligs,rnas: rrna, proteins: rps}}
    """.format_map({"pdbid":pdbid})
    return Response(_neoget(cypher))






@api_view(['GET', 'POST'])
def match_structs(request):
    params       = dict( request.GET )
    protstomatch = params['proteins'][0]
    targets      = protstomatch.split(',')
    targets      = [ *map(lambda x : f'\'{x}\'', targets) ]
    targets      = ",".join(targets)
    cypher       = """match (n:RibosomeStructure)-[]-(rp:Protein)
    with n, rp,[] as strnoms 
    unwind rp.nomenclature as unwound
    with collect(unwound) as unwound, n, [{targets}] as tgts
    where all(x in tgts where x in unwound)
    return n.rcsb_id""".format_map({"targets":targets})
    return Response(_neoget(cypher))

    
@api_view(['GET'])
def get_full_structure(request):
    params = dict(request.GET)
    pdbid  = str.upper(params['pdbid'][0])
    CYPHER_STRING="""
    match (rib:RibosomeStructure {{rcsb_id:'{}'}}) 
        unwind rib as rb
        optional match (l:Ligand)-[]-(rb)
        with collect(l.chemicalId) as ligs, rb
        optional match (rps:Protein)-[]-(rb)
        with ligs, rb, collect({{auth_asym_id:rps.auth_asym_id, nomenclature:rps.nomenclature, entity_poly_seq_one_letter_code: rps.entity_poly_seq_one_letter_code}}) as rps
        optional match (rnas:RNA)-[]-(rb)
        with ligs, rb, rps, collect({{auth_asym_id: rnas.auth_asym_id, nomenclature: rnas.nomenclature, entity_poly_seq_one_letter_code:rnas.entity_poly_seq_one_letter_code}}) as struct_rnas
        return {{
            struct : rb         ,
            ligands: ligs       ,
            rps    : rps        ,
            rnas   : struct_rnas
            }}
        """.format(pdbid)

    qres = _neoget(CYPHER_STRING)
    return Response(qres)
    
    

@api_view(['GET'])
def get_all_structs(request):
    CYPHER_STRING="""
    match (ribs:RibosomeStructure) 
        unwind ribs as rb
        optional match (l:Ligand)-[]-(rb)
        with collect(l.chemicalId) as ligs, rb
        optional match (rps:Protein)-[]-(rb)
        with ligs, rb, collect({{auth_asym_id:rps.auth_asym_id, nomenclature:rps.nomenclature, entity_poly_seq_one_letter_code: rps.entity_poly_seq_one_letter_code}}) as rps
        optional match (rnas:RNA)-[]-(rb)
        with ligs, rb, rps, collect({{auth_asym_id: rnas.auth_asym_id, nomenclature: rnas.nomenclature, entity_poly_seq_one_letter_code:rnas.entity_poly_seq_one_letter_code}}) as struct_rnas
        return {{
            struct : rb         ,
            ligands: ligs       ,
            rps    : rps        ,
            rnas   : struct_rnas
            }}
        """.format()

    qres = _neoget(CYPHER_STRING)
    return Response(qres)
#? ---------------------------PROTEINS

@api_view(['GET'])
def get_banclasses_metadata(request):
    params  = dict(request.GET)
    family  = str(params['family'][0]).lower()
    subunit = str(params['subunit'][0]).lower()


    if subunit == "ssu":
        fstring = 'toLower(n.class_id) contains "s" or toLower(n.class_id) contains "bthx" or toLower(n.class_id) contains "rack"' 
    elif subunit == "lsu": 
        fstring = 'toLower(n.class_id) contains "l"' 

    CYPHER_STRING="""
    match (n:ProteinClass)-[]-(rp:Protein)-[]-(s:RibosomeStructure) where  toLower(n.class_id) contains "{}"  and {} 
    unwind s.`src_organism_ids` as orgid
    with collect(distinct orgid) as allorgs, n as n, collect(s.rcsb_id) as structures, collect(distinct rp.pfam_comments) as comments
    return {{banClass: n.class_id, organisms:  allorgs, comments:comments, structs: structures }}""".format(family, fstring)

    return Response(_neoget(CYPHER_STRING))

@api_view(['GET'])
def list_nom_classes(request):
    CYPHER_STRING="""
    match (b:ProteinClass)-[]-(rp)-[]-(str:RibosomeStructure)
    with str, b, rp
    return {{
    nom_class: b.class_id,
    rps:collect({{
        organism_desc: rp.src_organism_names,
        organism_id  : rp.src_organism_ids,
        uniprot      : rp.uniprot_accession,
        parent       : str.rcsb_id,
        parent_reso  : str.resolution,
        strand_id    : rp.entity_poly_strand_id
        }}),
    presentIn:collect(str.rcsb_id)}}""".format_map({})
    return Response(_neoget(CYPHER_STRING))

@api_view(['GET','POST'])
def gmo_nom_class(request):
    params = dict(request.GET)
    ban    = str(params['banName'][0])
    CYPHER_STRING="""
    match (rib:RibosomeStructure)-[]-(n:Protein)-[]-(nc:ProteinClass{{class_id:"{ban}"}})

    return {{  
    parent_resolution                  : rib.resolution,
    parent_year                        : rib.citation_year,
    parent_method                      : rib.expMethod,
    parent_rcsb_id                     : n.parent_rcsb_id,
    pfam_accessions                    : n.pfam_accessions,
    pfam_comments                      : n.pfam_comments,
    pfam_descriptions                  : n.pfam_descriptions,
    src_organism_ids                   : n.src_organism_ids,
    src_organism_names                 : n.src_organism_names,
    uniprot_accession                  : n.uniprot_accession,
    rcsb_pdbx_description              : n.rcsb_pdbx_description,
    entity_poly_strand_id              : n.entity_poly_strand_id,
    entity_poly_seq_one_letter_code    : n.entity_poly_seq_one_letter_code,
    entity_poly_seq_one_letter_code_can: n.entity_poly_seq_one_letter_code_can,
    entity_poly_seq_length             : n.entity_poly_seq_length,
    entity_poly_polymer_type           : n.entity_poly_polymer_type,
    entity_poly_entity_type            : n.entity_poly_entity_type,
    surface_ratio                      : n.surface_ratio,
    nomenclature                       : n.nomenclature
        }}
    """.format_map({
        "ban":ban
    })

    return Response(_neoget(CYPHER_STRING))

@api_view(['GET']) 
def banclass_annotation(request):
    params   = dict(request.GET)
    banclass = str(params['banclass'][0])

    CYPHER_STRING=f"""
    
            match (n:ProteinClass{{class_id:"{banclass}"}})-[]-(rp:Protein) 
            with rp.rcsb_pdbx_description as dd 
            return  dd limit 6;
           
           """
    return Response(_neoget(CYPHER_STRING))

@api_view(['GET'])
def nomclass_visualize(request):

    params = dict(request.GET)
    ban    = str(params['ban'][0])

    CYPHER_STRING="""
    match (n:ProteinClass)-[]-(rp:Protein) where n.class_id ="{}" return  {{
    class: n.class_id,
    members: collect({{parent: rp.parent_rcsb_id, chain:rp.entity_poly_strand_id}}),
    comments:collect(distinct rp.pfam_comments)}} """.format(ban)
    
    return Response(_neoget(CYPHER_STRING))

@api_view(['GET'])
def proteins_number(request):
    CYPHER_STRING="""match (n:Protein) return count(n);"""
    return Response(_neoget(CYPHER_STRING))


#? ---------------------------RNA
@api_view(['GET']) 
def get_rnas_by_struct(request):
    CYPHER_STRING="""match (n:RibosomeStructure)-[]-(r:RNA) where toLower(r.rcsb_pdbx_description)  contains "mrna" or toLower(r.rcsb_pdbx_description) contains "trna"
        or toLower(r.rcsb_pdbx_description)  contains "m-rna" or toLower(r.rcsb_pdbx_description)  contains "t-rna"  or toLower(r.rcsb_pdbx_description)  contains "messenger" or toLower(r.rcsb_pdbx_description)  contains "transfer"
        return {{struct:n.rcsb_id, rnas:collect(r.rcsb_pdbx_description)}};""".format_map({})
    return Response(_neoget(CYPHER_STRING))


@api_view(['GET']) 
def get_rna_class(request):

    params        = dict(request.GET)
    rna_class     = str(params['rna_class'][0])


    CYPHER_STRING  = """
    match (c:RNAClass {{ class_id:"{}" }})-[]-(n)-[]-(rib:RibosomeStructure)
    return {{
        parent_year                         : rib.citation_year                      ,
        parent_resolution                   : rib.resolution                         ,
        parent_citation                     : rib.citation_title                     ,
        parent_method                       : rib.expMethod                          ,
        asym_ids                            : n.asym_ids                           ,
        auth_asym_id                        : n.auth_asym_id                       ,
        nomenclature                        : c.class_id                           ,
        parent_rcsb_id                      : n.parent_rcsb_id                     ,
        src_organism_names                  : n.src_organism_names                 ,
        host_organism_names                 : n.host_organism_names                ,
        src_organism_ids                    : n.src_organism_ids                   ,
        host_organism_ids                   : n.host_organism_ids                  ,
        rcsb_pdbx_description               : n.rcsb_pdbx_description              ,
        entity_poly_strand_id               : n.entity_poly_strand_id              ,
        entity_poly_seq_one_letter_code     : n.entity_poly_seq_one_letter_code    ,
        entity_poly_seq_one_letter_code_can : n.entity_poly_seq_one_letter_code_can,
        entity_poly_seq_length              : n.entity_poly_seq_length             ,
        entity_poly_polymer_type            : n.entity_poly_polymer_type           ,
        entity_poly_entity_type             : n.entity_poly_entity_type            ,
        ligand_like                         : n.ligand_like                        
    }}
    """.format(rna_class)

    return Response(_neoget(CYPHER_STRING))

@api_view(['GET'])
def anything(request):
    return Response("This a testing endpoint")

@api_view(['GET'])
def nomenclature(request):
    params        = dict(request.GET)
    if 'rcsb_id' in params and len( params['rcsb_id'] ) > 0:
        cypher_single ="""match (n:RibosomeStructure{{rcsb_id:"{}"}})-[]-(c) where c:Protein or c:RNA 
        return {{struct:n.rcsb_id, strand:c.auth_asym_id,nomenclature:c.nomenclature}}""".format(params['rcsb_id'][0].upper())
        maps=  _neoget(cypher_single)
        d={}
        for _ in maps:
            d.update({ _['auth_asym_id']:_['nomenclature'] })
        return Response(d)
    else:
        cypher_all ="""match (n:RibosomeStructure)-[]-(c) where c:Protein or c:RNA 
        return {{struct:n.rcsb_id, auth_asym_id:c.auth_asym_id,nomenclature:c.nomenclature}}""".format()
        all_maps = {

        }

        resp = _neoget(cypher_all)
        for _ in resp:
            if _['struct'] in all_maps:
                all_maps[_['struct']].update({_['auth_asym_id']:_['nomenclature']})
            else:
                all_maps.update(
                    {
                    _['struct']:{_['auth_asym_id']:_['nomenclature']}
                    }
                )

        return Response(all_maps)

def tax_ids(request):
    CYPHER_STRING="""
    match (r:RibosomeStructure) 
    unwind r.src_organism_ids as orgs
    return  collect(distinct orgs);
    """
    qres = _neoget(CYPHER_STRING)
    BY_STRUCT_CYPHER="""
    match (r:RibosomeStructure) 
    unwind r.src_organism_ids as orgs
    with distinct orgs as orgs
    match (s:RibosomeStructure) where orgs in s.`src_organism_ids`
    return {{organism: orgs, struct: s.rcsb_id}}  limit 5000
    """.format()
    by_struct = _neoget(BY_STRUCT_CYPHER)
    d         = {}
    for _ in by_struct:
        if _['organism'] not in d:
            d[_['organism']] = [_['struct']]
        else:
            d[_['organism']].append(_['struct'])
    return Response(d)

@api_view(['GET','POST'])
def custom_cypher(request):
    params        = dict(request.GET)
    CYPHER_STRING = params['cypher'][0]
    print("GOT STRING", CYPHER_STRING)
    k = _neoget(CYPHER_STRING)
    print("->>>>>>>>>>>>>",k)
    return Response()

#? ------------------------------ General 

# @api_view(['GET', 'POST'])
# def cif_chain(request):
#     params   = dict(request.GET)
#     structid = str(params['structid'])[0]
#     chainid  = str(params['chainid'])[0]

