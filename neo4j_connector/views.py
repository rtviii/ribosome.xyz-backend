from operator import sub
from rest_framework.decorators import api_view
from rest_framework.response import Response
from neo4j import  Result, GraphDatabase
import os
#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯
uri         =  os.environ[ 'NEO4J_URI' ]
authglobal  =  (os.environ[ 'NEO4J_USER' ],os.environ[ 'NEO4J_PASSWORD' ])
current_db  =  os.environ[ 'NEO4J_CURRENTDB' ]

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

# @api_view(['GET', 'POST'])
# def struct_to_ligands(request):
    
#     params = dict(request.GET)
#     struct_id = params['struct_id'][0]
#     cypher ="""
#     match (n:Ligand{{chemicalId:"{}"}}) return n;""".format(chemId)
#     return Response(_neoget(cypher))

@api_view(['GET', 'POST'])
def get_all_ligands(request):
    CYPHER_STRING="""
        match (l:Ligand)-[]-(r:RibosomeStructure)
        return {{
             ligand: l, 
         presentIn:collect({{
            _organismId          : r._organismId,                 
            rcsb_id              : r.rcsb_id,                 
            expMethod            : r.expMethod,                 
            resolution           : r.resolution,                 
            citation_title       : r.citation_title
            }})      
            }}
         
    """.format_map({})
    return Response(_neoget(CYPHER_STRING))


#? ---------------------------STRUCTS

@api_view(['GET']) 
def get_ligands_by_struct(request):
    CYPHER_STRING="""match (n:RibosomeStructure)-[]-(l:Ligand)
           return {{ title: n.citation_title, struct:n.rcsb_id, organism:n._organismName, taxid:n._organismId, 
           ligands:collect({{ chemid: l.chemicalId, name:l.chemicalName }})}}""".format_map({})
    return Response(_neoget(CYPHER_STRING))

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
    return {{structure: n, ligands: ligs,rnas: rrna, rps: rps}}
    """.format_map({"pdbid":pdbid})




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
#? ---------------------------PROTEINS

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
    match (n:RPClass)-[]-(rp:RibosomalProtein)-[]-(s:RibosomeStructure) where  toLower(n.class_id) contains "{}"  and {} 
    unwind s.`_organismId` as orgid
    with collect(distinct orgid) as allorgs, n as n, collect(s.rcsb_id) as structures, collect(distinct rp.pfam_comments) as comments
    return {{banClass: n.class_id, organisms:  allorgs, comments:comments, structs: structures }}""".format(family, fstring)

    return Response(_neoget(CYPHER_STRING))

@api_view(['GET'])
def list_nom_classes(request):
    CYPHER_STRING="""
    match (b:RPClass)-[]-(rp)-[]-(str:RibosomeStructure)
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
    match (rib:RibosomeStructure)-[]-(n:RibosomalProtein)-[]-(nc:RPClass{{class_id:"{ban}"}})

    return {{  parent_resolution                  : rib.resolution,
    parent_year                        : rib.citation_year,
    parent_method                      : rib.expMethod,
    parent_rcsb_id                      :n.parent_rcsb_id,
    pfam_accessions                     :n.pfam_accessions,
    pfam_comments                       :n.pfam_comments,
    pfam_descriptions                   :n.pfam_descriptions,
    rcsb_source_organism_id             :n.rcsb_source_organism_id,
    rcsb_source_organism_description    :n.rcsb_source_organism_description,
    uniprot_accession                   :n.uniprot_accession,
    rcsb_pdbx_description               :n.rcsb_pdbx_description,
    entity_poly_strand_id               :n.entity_poly_strand_id,
    entity_poly_seq_one_letter_code     :n.entity_poly_seq_one_letter_code,
    entity_poly_seq_one_letter_code_can :n.entity_poly_seq_one_letter_code_can,
    entity_poly_seq_length              :n.entity_poly_seq_length,
    entity_poly_polymer_type            :n.entity_poly_polymer_type,
    entity_poly_entity_type             :n.entity_poly_entity_type,
    surface_ratio                       :n.surface_ratio,
    nomenclature                        :n.nomenclature
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
    
            match (n:RPClass{{class_id:"{banclass}"}})-[]-(rp:RibosomalProtein) 
            with rp.rcsb_pdbx_description as dd 
            return  dd limit 6;
           
           """
    return Response(_neoget(CYPHER_STRING))

@api_view(['GET'])
def nomclass_visualize(request):

    params = dict(request.GET)
    ban    = str(params['ban'][0])

    CYPHER_STRING="""
    match (n:RPClass)-[]-(rp:RibosomalProtein) where n.class_id ="{}" return  {{
    class: n.class_id,
    members: collect({{parent: rp.parent_rcsb_id, chain:rp.entity_poly_strand_id}}),
    comments:collect(distinct rp.pfam_comments)}} """.format(ban)
    
    return Response(_neoget(CYPHER_STRING))


@api_view(['GET'])
def proteins_number(request):
    CYPHER_STRING="""match (n:RibosomalProtein) return count(n);"""
    return Response(_neoget(CYPHER_STRING))


#? ---------------------------RNA
@api_view(['GET']) 
def get_rnas_by_struct(request):
    CYPHER_STRING="""match (n:RibosomeStructure)-[]-(r:rRNA) where toLower(r.rcsb_pdbx_description)  contains "mrna" or toLower(r.rcsb_pdbx_description) contains "trna"
        or toLower(r.rcsb_pdbx_description)  contains "m-rna" or toLower(r.rcsb_pdbx_description)  contains "t-rna"  or toLower(r.rcsb_pdbx_description)  contains "messenger" or toLower(r.rcsb_pdbx_description)  contains "transfer"
        return {{struct:n.rcsb_id, rnas:collect(r.rcsb_pdbx_description)}};""".format_map({})
    return Response(_neoget(CYPHER_STRING))


@api_view(['GET']) 
def get_rna_class(request):
    params        = dict(request.GET)
    rna_class     = str(params['rna_class'][0])


    param2class = {
        '5'   : '5SrRNA',
        '5.8' : '5.8SrRNA',
        '12'  : '12SrRNA',
        '16'  : '16SrRNA',
        '21'  : '21SrRNA',
        '23'  : '23SrRNA',
        '25'  : '25SrRNA',
        '28'  : '28SrRNA',
        '35'  : '35SrRNA',
        'mrna': 'mRNA',
        'trna': 'tRNA',
    }





    CYPHER_STRING  = """match (c:RNAClass {{ class_id:"{}" }})-[]-(n)-[]-(rib:RibosomeStructure)
    return {{
    struct           : n.parent_rcsb_id,
    orgid            : n.rcsb_source_organism_id,
    description      : n.rcsb_pdbx_description,
    seq              : n.entity_poly_seq_one_letter_code,
    strand           : n.entity_poly_strand_id,
    parent_year      : rib.citation_year,
    parent_resolution: rib.resolution,
    parent_citation  : rib.citation_title,
    parent_method    : rib.expMethod,
    nomenclature     : c.class_id}}
    """.format(param2class[rna_class])


    return Response(_neoget(CYPHER_STRING))

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


    
    


@api_view(['GET','POST'])
def custom_cypher(request):
    params        = dict(request.GET)
    CYPHER_STRING = params['cypher']
    return Response(_neoget(CYPHER_STRING))