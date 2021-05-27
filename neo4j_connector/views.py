from operator import sub
from rest_framework.decorators import api_view
from rest_framework.response import Response
from neo4j import  Result, GraphDatabase
import os
#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯
uri         =  os.getenv( 'NEO4J_URI' )
authglobal  =  (os.getenv( 'NEO4J_USER' ),os.getenv( 'NEO4J_PASSWORD' ))
current_db  =  os.getenv( 'NEO4J_CURRENTDB' )
#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯

def _neoget(CYPHER_STRING:str):
    driver = GraphDatabase.driver(uri, auth= authglobal )
    def parametrized_query(tx, **kwargs):
        result:Result = tx.run(CYPHER_STRING, **kwargs)
        return result.value()

    with driver.session() as session:
        return session.read_transaction(parametrized_query)

#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯
@api_view(['GET'])
def anything(request):
    return Response("This a testing endpoint")


@api_view(['GET'])
def TEMP_classification_sample(request):
    CYPHER_STRING="""
    match (n:RibosomalProtein)-[]-(f:RibosomeStructure) where n.surface_ratio is not null 
    with distinct f as rb
        optional match (l:Ligand)-[]-(rb)
        with collect(l.chemicalId) as ligs, rb
        optional match (rps:RibosomalProtein)-[]-(rb)
        with ligs, rb, collect({{strands:rps.entity_poly_strand_id,surface_ratio:rps.surface_ratio, noms:rps.nomenclature}}) as rps
        optional match (rnas:rRNA)-[]-(rb)
        with ligs, rb, rps, collect(rnas.entity_poly_strand_id) as rnas
        return {{struct:rb, ligands: ligs, rps:rps, rnas:rnas}}""".format()
    qres = _neoget(CYPHER_STRING)
    return Response(qres)


@api_view(['GET'])
def get_all_structs(request):
    CYPHER_STRING="""
    
    match (ribs:RibosomeStructure) 
        unwind ribs as rb
        optional match (l:Ligand)-[]-(rb)
        with collect(l.chemicalId) as ligs, rb
        optional match (rps:RibosomalProtein)-[]-(rb)
        with ligs, rb, collect({{strands:rps.entity_poly_strand_id,surface_ratio:rps.surface_ratio, noms:rps.nomenclature}}) as rps
        optional match (rnas:rRNA)-[]-(rb)
        with ligs, rb, rps, collect(rnas.entity_poly_strand_id) as rnas
        return {{struct:rb, ligands: ligs, rps:rps, rnas:rnas}}
        """.format()

    qres = _neoget(CYPHER_STRING)
    return Response(qres)

@api_view(['GET'])
def get_surface_ratios(request):
    params = dict(request.GET)
    struct = params['pdbid'][0]
    cypher ="""match (r:RibosomeStructure{{rcsb_id:"{}"}})-[]-(p:RibosomalProtein)
    return {{strand:p.entity_poly_strand_id,nom:p.nomenclature,ratio:p.surface_ratio }}""".format(struct)
    return Response(_neoget(cypher))
    
@api_view(['GET', 'POST'])
def get_individual_ligand(request):
    params = dict(request.GET)
    chemId = params['chemId'][0]
    cypher ="""
    match (n:Ligand{{chemicalId:"{}"}}) return n;""".format(chemId)
    return Response(_neoget(cypher))
    

@api_view(['GET', 'POST'])
def match_structs(request):
    params       = dict( request.GET )
    protstomatch = params['proteins'][0]
    targets      = protstomatch.split(',')
    targets      = [ *map(lambda x : f'\'{x}\'', targets) ]
    targets      = ",".join(targets)
    cypher       = """match (n:RibosomeStructure)-[]-(rp:RibosomalProtein)
    with n, rp,[] as strnoms 
    unwind rp.nomenclature as unwound
    with collect(unwound) as unwound, n, [{targets}] as tgts
    where all(x in tgts where x in unwound)
    return n.rcsb_id""".format_map({"targets":targets})

    return Response(_neoget(cypher))

@api_view(['GET'])
def test_endpoint(request):
    return Response("This is testing endpoint. What did you expect?")

@api_view(['GET'])
def get_struct(request):
    params = dict(request.GET)
    pdbid  = str.upper(params['pdbid'][0])
    cypher = """
    match (n:RibosomeStructure{{rcsb_id:"{pdbid}"}})
    optional match (rr:rRNA)-[]-(n)
    with n, collect(rr) as rrna
    optional match (rp:RibosomalProtein)-[]-(n)
    with n, rrna,  collect(rp) as rps
    optional match (l:Ligand)-[]-(n)
    with n, rrna, rps, collect(l) as ligs
    return {{structure: n, ligands:ligs,rnas:rrna,rps:rps}}""".format_map({"pdbid":pdbid})


    return Response(_neoget(cypher))


@api_view(['GET'])
def get_homologs(request):
    params        = dict(request.GET)
    ban           = str(params['banName'][0])

    CYPHER_STRING = """
    match (s:RibosomalProtein)-[]-(q:RibosomeStructure) where 
    "{ban}" in s.nomenclature  
    return {{protein:s, subchain_of: q.`_PDBId`}}""".format_map({
            "ban":ban
        })
    return Response(_neoget(CYPHER_STRING))



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
    match (n:NomenclatureClass)-[]-(rp:RibosomalProtein)-[]-(s:RibosomeStructure) where  toLower(n.class_id) contains "{}"  and {} 
    unwind s.`_organismId` as orgid
    with collect(distinct orgid) as allorgs, n as n, collect(s.rcsb_id) as structures, collect(distinct rp.pfam_comments) as comments
    return {{banClass: n.class_id, organisms:  allorgs, comments:comments, structs: structures }}""".format(family, fstring)


    return Response(_neoget(CYPHER_STRING))

@api_view(['GET'])
def list_nom_classes(request):
    CYPHER_STRING="""
    match (b:NomenclatureClass)-[]-(rp)-[]-(str:RibosomeStructure)
    with str, b, rp
    return {{
    nom_class: b.class_id,
    rps:collect({{
        organism_desc: rp.rcsb_source_organism_description,
        organism_id  : rp.rcsb_source_organism_id,
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
    match (q:RibosomeStructure)-[]-(n:RibosomalProtein)-[]-(nc:NomenclatureClass{{class_id:"{ban}"}})
    return collect(n);
    """.format_map({
        "ban":ban
    })

    return Response(_neoget(CYPHER_STRING))


@api_view(['GET', 'POST'])
def get_all_ligands(request):
    CYPHER_STRING="""

        match (l:Ligand)-[]-(r:RibosomeStructure)
        return {{
             ligand: l, 
         presentIn:collect({{
            _organismId          : r._organismId,                 
            _organismName        : r._organismName,                 
            rcsb_id              : r.rcsb_id,                 
            expMethod            : r.expMethod,                 
            resolution           : r.resolution,                 
            cryoem_exp_resolution: r.cryoem_exp_resolution,                 
            citation_title       : r.citation_title,                 
            pdbx_keywords_text   : r.pdbx_keywords_text}})      }}


         
    """.format_map({})
    return Response(_neoget(CYPHER_STRING))


@api_view(['GET']) 
def get_rna_class(request):
    params        = dict(request.GET)
    rna_class     = str(params['rna_class'][0])
    CYPHER_STRING = """
    match (n:rRNA) where toLower(n.rcsb_pdbx_description) contains '{}'
    return {{
    struct     : n.parent_rcsb_id,
    orgid      : n.rcsb_source_organism_id,
    description: n.rcsb_pdbx_description,
    seq        : n.entity_poly_seq_one_letter_code,
        strand: n.entity_poly_strand_id 
    
    }}
    """.format(rna_class)

    if rna_class == '5':
        CYPHER_STRING = """
        match (n:rRNA) where toLower(n.rcsb_pdbx_description) contains '5' and  
        not ( toLower(n.rcsb_pdbx_description) contains '5.8'
        or toLower(n.rcsb_pdbx_description) contains '4.55'
        or toLower(n.rcsb_pdbx_description) contains '25'
        or toLower(n.rcsb_pdbx_description) contains '35'
        )
        return {{
        struct     : n.parent_rcsb_id,
        orgid      : n.rcsb_source_organism_id,
        description: n.rcsb_pdbx_description,
        seq        : n.entity_poly_seq_one_letter_code,
        strand     : n.entity_poly_strand_id
        }}
        """.format(rna_class)
    if rna_class == 'trna':
        CYPHER_STRING = """
        match (n:rRNA) where toLower(n.rcsb_pdbx_description) contains 'trna' or  toLower(n.rcsb_pdbx_description) contains 'transport'
        return {{
        struct     : n.parent_rcsb_id,
        orgid      : n.rcsb_source_organism_id,
        description: n.rcsb_pdbx_description,
        seq        : n.entity_poly_seq_one_letter_code,
        strand: n.entity_poly_strand_id
        }}
        """.format(rna_class)
    if rna_class == 'mrna':
        CYPHER_STRING = """
        match (n:rRNA) where toLower(n.rcsb_pdbx_description) contains 'mrna' or  toLower(n.rcsb_pdbx_description) contains 'messenger'
        return {{
        struct     : n.parent_rcsb_id,
        orgid      : n.rcsb_source_organism_id,
        description: n.rcsb_pdbx_description,
        seq        : n.entity_poly_seq_one_letter_code,
        strand: n.entity_poly_strand_id
        }}
        """.format(rna_class)
    if rna_class == 'other':

        CYPHER_STRING = """match (n:rRNA) where not
        (toLower(n.rcsb_pdbx_description) contains   '5'
        or toLower(n.rcsb_pdbx_description) contains '5.8'
        or toLower(n.rcsb_pdbx_description) contains '12'
        or toLower(n.rcsb_pdbx_description) contains '21'
        or toLower(n.rcsb_pdbx_description) contains '16'
        or toLower(n.rcsb_pdbx_description) contains '18'
        or toLower(n.rcsb_pdbx_description) contains '23'
        or toLower(n.rcsb_pdbx_description) contains '26'
        or toLower(n.rcsb_pdbx_description) contains '28'
        or toLower(n.rcsb_pdbx_description) contains 'trna'
        or toLower(n.rcsb_pdbx_description) contains 'transport'
        or toLower(n.rcsb_pdbx_description) contains 'mrna'
        or toLower(n.rcsb_pdbx_description) contains 'messenger'
        )

        return {
        struct     : n.parent_rcsb_id,
        orgid      : n.rcsb_source_organism_id,
        description: n.rcsb_pdbx_description,
        seq        : n.entity_poly_seq_one_letter_code,
        strand: n.entity_poly_strand_id
        }

        """

    return Response(_neoget(CYPHER_STRING))


@api_view(['GET','POST'])
def custom_cypher(request):
    params        = dict(request.GET)
    CYPHER_STRING = params['cypher']
    return Response(_neoget(CYPHER_STRING))


@api_view(['GET']) 
def get_rnas_by_struct(request):
    CYPHER_STRING="""match (n:RibosomeStructure)-[]-(r:rRNA) where toLower(r.rcsb_pdbx_description)  contains "mrna" or toLower(r.rcsb_pdbx_description) contains "trna"
        or toLower(r.rcsb_pdbx_description)  contains "m-rna" or toLower(r.rcsb_pdbx_description)  contains "t-rna"  or toLower(r.rcsb_pdbx_description)  contains "messenger" or toLower(r.rcsb_pdbx_description)  contains "transfer"
        return {{struct:n.rcsb_id, rnas:collect(r.rcsb_pdbx_description)}};""".format_map({})
    return Response(_neoget(CYPHER_STRING))

@api_view(['GET']) 
def get_ligands_by_struct(request):
    CYPHER_STRING="""match (n:RibosomeStructure)-[]-(l:Ligand)
           return {{ title: n.citation_title, struct:n.rcsb_id, organism:n._organismName, taxid:n._organismId, 
           ligands:collect({{ chemid: l.chemicalId, name:l.chemicalName }})}}""".format_map({})
    return Response(_neoget(CYPHER_STRING))


@api_view(['GET'])
def nomclass_visualize(request):

    params = dict(request.GET)
    ban    = str(params['ban'][0])

    CYPHER_STRING="""
    match (n:NomenclatureClass)-[]-(rp:RibosomalProtein) where n.class_id ="{}" return  {{
    class: n.class_id,
    members: collect({{parent: rp.parent_rcsb_id, chain:rp.entity_poly_strand_id}}),
    comments:collect(distinct rp.pfam_comments)}} """.format(ban)
    
    return Response(_neoget(CYPHER_STRING))