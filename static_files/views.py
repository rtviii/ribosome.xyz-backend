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

from rxz_backend.settings import PROJECT_PATH


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

@api_view(['GET','POST'])
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

@api_view(['GET','POST'])
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

@api_view(['GET','POST'])
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

@api_view(['GET','POST'])
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

@api_view(['GET','POST'])
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

@api_view(['GET', 'POST'])
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

@api_view(['GET', 'POST'])
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

@api_view(['GET', 'POST'])
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

@api_view(['GET', 'POST'])
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

@api_view(['GET', 'POST'])
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

