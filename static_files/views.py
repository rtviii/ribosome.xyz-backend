from cgi import parse_multipart
from genericpath import isfile
from multiprocessing import ProcessError
from os import error
import sys
from django.template import response
from dotenv import load_dotenv
from io import StringIO, BytesIO
from neo4j import GraphDatabase
from rest_framework.decorators import api_view
from rest_framework.response import Response
import os
import tempfile
from rxz_backend.settings import STATIC_ROOT
import json
from django.http import FileResponse, HttpResponse
from wsgiref.util import FileWrapper
import zipfile
from neo4j import  Result, GraphDatabase
import subprocess
from ribetl.ciftools import transpose_ligand


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



@api_view(['GET','POST'])
def pairwise_align(request):
    params = dict(request.GET)

    print("-------------------+------------------")
    struct1 = params['struct1'][0].upper()
    struct2 = params['struct2'][0].upper()

    strand1 = params['strand1'][0]
    strand2 = params['strand2'][0]

    name1   = "{}_STRAND_{}.cif".format(struct1,strand1)
    name2   = "{}_STRAND_{}.cif".format(struct2,strand2)
    print(f"Attempting to align \033[95m {name1}\033[0m with \033[95m{name2}\033[0m.")
    handle1 = os.path.join(STATIC_ROOT, struct1, "CHAINS", name1)
    handle2 = os.path.join(STATIC_ROOT, struct2, "CHAINS", name2)

    for x in [handle1,handle2]:
        if not os.path.isfile(x):
            raise FileNotFoundError(f"File {x} is not found in {STATIC_ROOT}")
    
    protein_alignment_script = os.environ['PROTEIN_ALIGNMENT_SCRIPT']
    print("Using alignment script ", protein_alignment_script)

    # subprocess.call(f'{protein_alignment_script} {handle1} {handle2} {struct1+"_"+struct2} {struct2+"_"+strand2}')
    try:
        subprocess.call([protein_alignment_script, handle1 ,handle2 ,struct1+"_"+strand1, struct2+"_"+strand2])
    except:
        raise ProcessError

    # os.system(f'{protein_alignment_script} {handle1} {handle2} {struct1+"_"+struct2} {struct2+"_"+strand2}')

    alignedfile=os.environ["TEMP_CHAIN"]
    print("Produced alignment file:", alignedfile)

    try:
        doc = open(alignedfile, 'rb')
    except: 
        print(f"Could not find {alignedfile}. Exited")
        return Response(-1)

    response = HttpResponse(FileWrapper(doc), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}-{}_{}-{}.cif"'.format(
    struct1,strand1,struct2,strand2)

    return response

def fetch_strand(structid:str,strandid:str)->FileWrapper:
    filename   = "{}_STRAND_{}.cif".format(structid.upper(),strandid)
    filehandle = os.path.join(STATIC_ROOT, structid.upper(),'CHAINS', filename)

    # try: 
    doc = open(filehandle)
    return doc


@api_view(['GET','POST'])
def get_static_catalogue(request):
    file_handle = os.path.join(os.getenv('STATIC_ROOT'),'static_files_catalogue.json')
    with open(file_handle) as infile:
        catalogue = json.load(infile)
        print(catalogue)
    return Response(catalogue)


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
    params = dict(request.GET)
    src_struct = params['src_struct' ][0].upper()
    chemid     = params['chemid'     ][0].upper()
    tgt_struct = params['tgt_struct' ][0].upper()

    print("Attempting to render ligand {} from {}(orig) in {}.".format(chemid,src_struct,tgt_struct))
    prediction_filename   = "PREDICTION_{}_{}_{}.json".format(chemid,src_struct,tgt_struct)
    filehandle = os.path.join(STATIC_ROOT, tgt_struct, prediction_filename)

    # print(">>>>>>>>>>>>>>>>>>>\033[93m Attempting to open \033[0m", filehandle)
    # print("got params", filehandle)

    #* Transpose Ligand Script

    transpose_ligand.init_transpose_ligand(src_struct,tgt_struct,chemid)


    try:
        with open(filehandle) as infile:
            data = json.load(infile)
            return Response(data)
    except error: 
        return Response(-1)


@api_view(['GET','POST'])
def get_ligand_nbhd(request):
    params = dict(request.GET)
    struct = params['struct'][0].upper()
    chemid = params['chemid'][0].upper()

    filename   = "LIGAND_{}.json".format(chemid)
    filehandle = os.path.join(STATIC_ROOT, struct, filename)
    print(">>>>>>>>>>>>>>>>>>>\033[93m Attempting to open \033[0m", filehandle)


    print("got params", filehandle)
    try:
        with open(filehandle) as infile:

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

    filename   = "{}.cif".format(struct_id)
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


    CYPHER = """match (n:RibosomeStructure)-[]-(r:RibosomalProtein)-[]-(b:NomenclatureClass)
    where n.rcsb_id ="{}" and b.class_id = "{}"
    return {{ struct: n.rcsb_id, strand: r.entity_poly_strand_id }}""".format(struct,classid)
    print("C----->" ,CYPHER)
    chains = _neoget(CYPHER)
    
    if len( chains ) < 1 :
        return Response(NotFoundErr)
    strand = chains[0]['strand']
    filename   = "{}_STRAND_{}.cif".format(struct,strand)
    filehandle = os.path.join(STATIC_ROOT, struct,'CHAINS', filename)


    print(filehandle)
    try:
        doc = open(filehandle, 'rb')
    except: 
        return Response("File not found")

    response = HttpResponse(FileWrapper(doc), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    return response

    # return Response(chains[0]['strand'])
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

