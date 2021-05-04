from os import error
import sys
from dotenv import load_dotenv
from rest_framework.decorators import api_view
from rest_framework.response import Response
import os
from rxz_backend.settings import STATIC_ROOT
import json
from django.http import HttpResponse
from wsgiref.util import FileWrapper
import zipfile
import subprocess


@api_view(['GET','POST'])
def pairwise_align(request):
    params = dict(request.GET)

    struct1 = params['struct1'][0].upper()
    struct2 = params['struct2'][0].upper()

    strand1 = params['strand1'][0]
    strand2 = params['strand2'][0]

    name1   = "{}_STRAND_{}.cif".format(struct1,strand1)
    name2   = "{}_STRAND_{}.cif".format(struct2,strand2)

    handle1 = os.path.join(STATIC_ROOT, struct1, "CHAINS", name1)
    handle2 = os.path.join(STATIC_ROOT, struct2, "CHAINS", name2)
    
    protein_alignment_script = os.getenv('PROTEIN_ALIGNMENT_SCRIPT')
    # subprocess.call(f'{protein_alignment_script} {handle1} {handle2} {struct1+"_"+struct2} {struct2+"_"+strand2}')
    subprocess.call([protein_alignment_script, handle1 ,handle2 ,struct1+"_"+strand1, struct2+"_"+strand2])

    # os.system(f'{protein_alignment_script} {handle1} {handle2} {struct1+"_"+struct2} {struct2+"_"+strand2}')

    alignedfile=os.getenv("TEMP_CHAIN")

    print("Alignment script:", protein_alignment_script)
    print("Alignment file:", alignedfile)
    try:
        doc = open(alignedfile, 'rb')
    except: 
        print(f"Could find {alignedfile}. Exited")
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
    params = dict(request.GET)
    structid = params['structid'][0].upper()
    chemid = params['chemid'][0].upper()

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
def get_ligand_nbhd(request):
    params = dict(request.GET)
    struct = params['struct'][0].upper()
    chemid = params['chemid'][0].upper()

    filename   = "LIGAND_{}.json".format(chemid)
    filehandle = os.path.join(STATIC_ROOT, struct, filename)


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
    params     =  request.GET.dict()
    file_names = [
        ]

    for item in params.items():
        file_names.append(os.path.join(STATIC_ROOT,item[0],'CHAINS',"{}_STRAND_{}.cif".format(item[0],item[1])))
    path         =  os.path.join(STATIC_ROOT, 'zipfile.zip')
    compression  =  zipfile.ZIP_DEFLATED

    # create the zip file first parameter path/name, second mode
    zf = zipfile.ZipFile(path, mode="w")
    try:
        for i,file_name in enumerate( file_names ):
            print(file_name)
            # Add file to the zip file
            # first parameter file to zip, second filename in zip
            zf.write( file_name, "{}_{}.cif".format(file_name,i), compress_type=compression)

    except FileNotFoundError:
        print("An error occurred")
    finally:
        # Don't forget to close the file!
        zf.close()

       
    response = HttpResponse(zf, content_type='application/force-download')
    response['Content-Disposition'] = f'attachment; filename="%s"' % 'compressed.zip'
    return response