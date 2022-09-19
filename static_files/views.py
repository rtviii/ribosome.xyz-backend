from os import error
from neo4j import GraphDatabase
from rest_framework.decorators import api_view
from rest_framework.response import Response
import os
import json
from django.http import FileResponse, HttpResponse
from wsgiref.util import FileWrapper
import zipfile
from neo4j import  Result, GraphDatabase
from ribetl.ciftools import transpose_ligand
from ribetl.ciftools.bsite_mixed import BindingSite

from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

from rest_framework.serializers import Serializer

uri         =  os.environ.get( 'NEO4J_URI'                                           )
authglobal  = (os.environ.get( 'NEO4J_USER'      ),os.environ.get( 'NEO4J_PASSWORD' ))
current_db  =  os.environ.get( 'NEO4J_CURRENTDB'                                     )
STATIC_ROOT =  os.environ.get( "STATIC_ROOT"                                         )
#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯

def _neoget(CYPHER_STRING:str):
    driver = GraphDatabase.driver(uri, auth= authglobal )
    def parametrized_query(tx, **kwargs):
        result:Result = tx.run(CYPHER_STRING, **kwargs)
        return result.value()

    with driver.session() as session:
        return session.read_transaction(parametrized_query)

r1star        = openapi.Parameter('r1star       ', openapi.IN_QUERY, description="Id of the FIRST residue of the FIRST of two given chains.", type=openapi.TYPE_STRING)
r1end         = openapi.Parameter('r1end        ', openapi.IN_QUERY, description="Id of the LAST residue of the FIRST of two given chains.", type=openapi.TYPE_STRING)
r2start       = openapi.Parameter('r2start      ', openapi.IN_QUERY, description="Id of the FIRST residue of the SECOND of two given chains.", type=openapi.TYPE_STRING)
r2end         = openapi.Parameter('r2end        ', openapi.IN_QUERY, description="Id of the LAST residue of the SECOND of two given chains.", type=openapi.TYPE_STRING)
struct1       = openapi.Parameter('struct1      ', openapi.IN_QUERY, description="Parent of the FIRST of the two aligned chains. ", type=openapi.TYPE_STRING)
struct2       = openapi.Parameter('struct2      ', openapi.IN_QUERY, description="Parent of the SECOND of the two aligned chains. ", type=openapi.TYPE_STRING)
auth_asym_id1 = openapi.Parameter('auth_asym_id1', openapi.IN_QUERY, description="auth_asym_id of the FIRST of the two aligned chains(according to RCSB)", type=openapi.TYPE_STRING)
auth_asym_id2 = openapi.Parameter('auth_asym_id2', openapi.IN_QUERY, description="auth_asym_id of the SECOND of the two aligned chains(according to RCSB)", type=openapi.TYPE_STRING)
@swagger_auto_schema(method='get',operation_description="Align substrings of two given chains in space and return the result as a x/mmcif file. The backbone of Superimposition tool. ", query_serializer=Serializer, manual_parameters=[
r1star,
r1end,
r2start,
r2end,
struct1,
struct2,
auth_asym_id1,
auth_asym_id2
])
@api_view(['GET',])
def ranged_align(request):
    params = dict(request.GET)
    print("-------------------+------------------")
    print("GOT PARAMS", params)
    print("-------------------+------------------")

    r1start = int(params['r1start'][0])
    r1end   = int(params['r1end'][0])

    r2start = int(params['r2start'][0])
    r2end   = int(params['r2end'][0])

    struct1       = params['struct1'][0].upper()
    struct2       = params['struct2'][0].upper()
    auth_asym_id1 = params['auth_asym_id1'][0]
    auth_asym_id2 = params['auth_asym_id2'][0]

    RANGED_ALIGNMENT_SCRIPT= os.path.join(str( PROJECT_PATH ), 'static_files','ranged_align.py')
    os.system("python3 {} {} {} {} {} {}-{} {}-{}".format(RANGED_ALIGNMENT_SCRIPT,struct1,struct2, auth_asym_id1, auth_asym_id2, r1start,r1end, r2start,r2end))

    alignedfile = os.environ["TEMP_CHAIN"]
    try:
        doc = open(alignedfile, 'rb')
    except: 
        print(f"Could not find {alignedfile}. Exited")
        return Response(-1)

    response = HttpResponse(FileWrapper(doc), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}-{}_{}-{}.cif"'.format(struct1,auth_asym_id1,struct2,auth_asym_id2)

    return response

get_chain_chainid = openapi.Parameter('chainid', openapi.IN_QUERY, description="auth_asym_id of the given chain(according to RCSB)", type=openapi.TYPE_STRING)
get_chain_structid = openapi.Parameter('structid', openapi.IN_QUERY, description="4-letter code of the parent structure.", type=openapi.TYPE_STRING)
@swagger_auto_schema(method='get',operation_description="Download a separate chain in a given structure as a .pdb file.", query_serializer=Serializer, manual_parameters=[
    get_chain_chainid,
    get_chain_structid
])
@api_view(['GET',])
def get_chain(request):
    params = dict(request.GET)
    chainid = params['chainid'][0]
    structid = str.upper(params['structid'][0])
    filename = "{}_subchain_{}.pdb".format(structid, chainid)

    file_handle = os.path.join(STATIC_ROOT,structid, filename)

    document = open(file_handle, 'rb')
    response = HttpResponse(FileWrapper(document), content_type='chemical/x-pdb')
    response['Content-Disposition'] = 'attachment; filename="{}_subchain_{}.pdb"'.format(structid, chainid)
    return response

@swagger_auto_schema(methods=[ 'get' ], auto_schema=None)  
@api_view(['GET',])
def download_ligand_nbhd(request):
    params   = dict(request.GET)
    structid = params['structid'][0].upper()
    chemid   = params['chemid'][0].upper()

    filename   = "LIGAND_{}.json".format(chemid)
    filehandle = os.path.join(STATIC_ROOT, structid, filename)

    try:
        doc = open(filehandle, 'rb')
    except: 
        print(f"Could find {filehandle}.Exited")
        return Response(-1)

    response = HttpResponse(FileWrapper(doc), content_type='application/json')
    response['Content-Disposition'] = 'attachment; filename="{}_LIGAND_{}.json"'.format(structid, chemid)
    return response

ligand_pred_src        = openapi.Parameter('src_struct', openapi.IN_QUERY, description="4-letter code of the SOURCE structure. Ex. '5AFI'. ", type=openapi.TYPE_STRING)
ligand_pred_tgt        = openapi.Parameter('tgt_struct', openapi.IN_QUERY, description="4-letter code of the TARGET structure. Ex. '7K00'. ", type=openapi.TYPE_STRING)
ligand_pred_liglike_id = openapi.Parameter('ligandlike_id', openapi.IN_QUERY, description="Chemical id of the ligand. Ex. \"PAR\" for Paromomycin", type=openapi.TYPE_STRING)
ligand_pred_is_polymer = openapi.Parameter('is_polymer', openapi.IN_QUERY, description="Whether the sought ligand is a polymer ( like an transcription factor or an mRNA) or a simple ligand. 'true' | 'false'", type=openapi.TYPE_STRING)
@swagger_auto_schema(method='get',operation_description="Download a prediction of given ligand's binding site in a target structure given an extant ligand and its parent structure.", query_serializer=Serializer, manual_parameters=[
ligand_pred_src        ,   
ligand_pred_tgt,
ligand_pred_liglike_id ,   
ligand_pred_is_polymer 
])
@api_view(['GET',])
def ligand_prediction(request):

    params     = dict(request.GET)
    # it is either a chemical id (if is_polymer == False) 
    # or a entity_poly_strand_id in the case of a ligand-like polymer (is_polymer  ==True)
    ligandlike_id = params['ligandlike_id'][0]
    src_struct    = params['src_struct' ][0].upper()
    tgt_struct    = params['tgt_struct' ][0].upper()
    is_polymer    = str( params['is_polymer' ][0] )
    print("Attempting to render  {} from {}(orig) in {}.".format(ligandlike_id,src_struct,tgt_struct))
    prediction_filename = "PREDICTION_{}_{}_{}.json".format     (ligandlike_id,src_struct ,tgt_struct          )
    filehandle          = os.path  .join(STATIC_ROOT, tgt_struct, prediction_filename)

    orig_binding_site_handle = os.path.join(STATIC_ROOT, src_struct, "{}_{}.json".format("POLYMER" if is_polymer.lower() == "true" else "LIGAND", ligandlike_id))
    target_json_handle       = os.path.join(STATIC_ROOT, tgt_struct, "{}.json".format(tgt_struct))
    with open(orig_binding_site_handle, 'rb') as infile:
        bsite = BindingSite(json.load(infile))

    with open(target_json_handle, 'rb') as target_handle:
        target_handle   = json.load(target_handle)

    #* Transpose Ligand Script
    prediction = transpose_ligand.init_transpose_ligand(src_struct,tgt_struct, target_handle, bsite)
    return Response(prediction)

ligand_nbhd_src        = openapi.Parameter('src_struct', openapi.IN_QUERY, description="4-letter code of the SOURCE structure. Ex. '5AFI'. ", type=openapi.TYPE_STRING)
ligand_nbhd_liglike_id = openapi.Parameter('ligandlike_id', openapi.IN_QUERY, description="Chemical id of the ligand. Ex. \"PAR\" for Paromomycin", type=openapi.TYPE_STRING)
ligand_nbhd_is_polymer = openapi.Parameter('is_polymer', openapi.IN_QUERY, description="Whether the sought ligand is a polymer ( like an transcription factor or an mRNA) or a simple ligand. 'true' | 'false'", type=openapi.TYPE_STRING)
@swagger_auto_schema(method='get',operation_description="Download a prediction of given ligand's binding site in a target structure given an extant ligand and its parent structure.", query_serializer=Serializer, 
                     manual_parameters=[
ligand_nbhd_src        ,
ligand_nbhd_liglike_id ,
ligand_nbhd_is_polymer 
])
@api_view(['GET',])
def get_ligand_nbhd(request):
    params        = dict(request.GET)
    print("----------------PARAMS ")
    print(params)
    src_struct    = params['src_struct'][0].upper()
    ligandlike_id = params['ligandlike_id'][0]
    is_polymer    = str( params['is_polymer'][0] )
    filehandle    = os.path .join(STATIC_ROOT, src_struct, "{}_{}.json".format("POLYMER" if is_polymer.lower() == 'true' else "LIGAND", ligandlike_id))

    print(f"Returning data from file {filehandle}")
    try:
        with open(filehandle, 'rb') as infile:
            data = json.load(infile)
            return Response(data)
    except error: 
        print("errored out", error)
        
        return Response(-1)




cif_chain_chainid        = openapi.Parameter('chainid', openapi.IN_QUERY, description="Auth_asym_id of a given chain according to RCSB.", type=openapi.TYPE_STRING)
cif_chain_structid        = openapi.Parameter('structid', openapi.IN_QUERY, description="4-letter code of the structure. Ex. '5AFI'. ", type=openapi.TYPE_STRING)
@swagger_auto_schema(method='get',operation_description="Download a given .mmcif chain.", query_serializer=Serializer, 
                     manual_parameters=[
                        cif_chain_chainid,
                        cif_chain_structid
])
@api_view(['GET', ])
def cif_chain(request):
    params     = dict(request.GET)
    chainid    = params['chainid'][0].upper()
    struct     = params['structid'][0].upper()

    filename   = "{}_STRAND_{}.cif".format(struct,chainid)
    filehandle = os.path.join(STATIC_ROOT, struct,'CHAINS', filename)

    try:
        doc = open(filehandle, 'rb')
    except: 
        return Response("File not found")

    response = HttpResponse(FileWrapper(doc), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)

    return response

@api_view(['GET'])
def download_structure(request):
    params     = dict(request.GET)
    struct_id    = params['struct_id'][0].upper()

    filename   = "{}_modified.cif".format(struct_id)
    filehandle = os.path.join(STATIC_ROOT, struct_id, filename)

    try:
        doc = open(filehandle, 'rb')
    except: 
        return Response("File not found")
    response = HttpResponse(FileWrapper(doc), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    return response

cif_chain_class_classid        = openapi.Parameter('classid', openapi.IN_QUERY, description="Subchain class of a chain to lookup: ex. uL4 or 35SrRNA.", type=openapi.TYPE_STRING)
cif_chain_class_struct        = openapi.Parameter('struct', openapi.IN_QUERY, description="4-letter code of the structure. Ex. '5AFI'. ", type=openapi.TYPE_STRING)
@swagger_auto_schema(method='get',operation_description="Attempt to find a chain of a given class and download the .mmcif file for it.", query_serializer=Serializer, 
                     manual_parameters=[
                        cif_chain_chainid,
                        cif_chain_structid
])
@api_view(['GET', ])
def cif_chain_by_class(request):
    params     = dict(request.GET)
    classid    = params['classid'][0]
    struct     = params['struct'][0].upper()

    CYPHER = """match (n:RibosomeStructure)-[]-(r)-[]-(b) where n.rcsb_id ="{}" and b.class_id = "{}"
    return {{ struct: n.rcsb_id, auth_asym_id: r.auth_asym_id }}""".format(struct,classid)

    chains = _neoget(CYPHER)
    if len( chains ) < 1 :
        return Response("Not found")
    auth_asym_id     = chains[0]['auth_asym_id']
    filename   = "{}_STRAND_{}.cif".format(struct,auth_asym_id)
    filehandle = os.path.join(STATIC_ROOT, struct.upper(), "CHAINS", filename)

    try:
        doc = open(filehandle, 'rb')
    except: 
        return Response("File not found")

    response = HttpResponse(FileWrapper(doc), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    return response


@swagger_auto_schema(methods=[ 'get' ], auto_schema=None)  
@api_view(['GET', ])
def tunnel(request):
    params     = dict(request.GET)
    struct     = params['struct'][0].upper()
    filetype   = params['filetype'][0]

    cetrline_filehandle  =  os.path.join(STATIC_ROOT, struct,'TUNNEL', 'csv', 'centerline.csv')
    report_filehandle    =  os.path.join(STATIC_ROOT, struct,'TUNNEL', f"{struct}_TUNNEL_REPORT.json")


    if filetype== 'report':
        print("GOT REQUEST FOR REPORT with params", params)
        try:
            doc = open(report_filehandle, 'rb')
        except: 
            return Response("File not found")

        response = HttpResponse(FileWrapper(doc), content_type='application/json')
        response['Content-Disposition'] = 'attachment; filename="{}"'.format(f"{struct}_tunnel_report.json")

        return response
    elif filetype == 'centerline':
        try:
            doc = open(cetrline_filehandle, 'rb')
        except: 
            return Response("File not found")

        response = HttpResponse(FileWrapper(doc), content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename="{}"'.format(f"{struct}_tunnel_centerline.csv")

        return response

@swagger_auto_schema(methods=[ 'get' ], auto_schema=None)  
@api_view(['GET', ])
def downloadArchive(request):

    params  = dict(request.GET)
    print("received dict " , params)

    if 'rna' in params:
        rnas    = params['rna']
        rnas  =[*map(lambda x: x.split('.'),rnas)]
        rnas  =[*map(lambda x:  os.path.join(STATIC_ROOT,x[0].upper(),'CHAINS',"{}_STRAND_{}.cif".format(x[0].upper(),x[1])),rnas)]
        print("dict:", rnas)
    else:
        rnas = []

    if 'structs' in params:
        structs = params['structs']
        structs = [*map(lambda x:  os.path.join(STATIC_ROOT,x.upper(),"{}.cif".format(x.upper())),structs)]
        print("structs", structs)
    else:
        structs = []

    file_names = [*structs,*rnas]
    zip_subdir = 'ribosome_xyz_archive'
    zf         = zipfile.ZipFile('temp_zip.zip', "w")

    for fpath in file_names:
        fdir, fname = os.path.split(fpath)
        zip_path = os.path.join(zip_subdir, fname)
        try:
            zf.write(fpath, zip_path)
            print(f"Failed to find file { fpath }")
        except:
            continue

    zf.close()
    r = open('temp_zip.zip','rb')
    os.remove('temp_zip.zip')
    return FileResponse(r)

