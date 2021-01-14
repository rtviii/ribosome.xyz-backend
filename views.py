from os import error
from rest_framework.decorators import api_view
from rest_framework.response import Response
import os, sys
from rxz_backend.settings import STATIC_ROOT
import json
from django.http import HttpResponse
from wsgiref.util import FileWrapper


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

    os.system('/home/rt/backend/static_files/align.sh {} {} {} {}'.format(
    handle1,handle2,
    struct1+"_"+strand1,
    struct2+"_"+strand2
    ))


    alignedfile='/home/rt/backend/static/TEMP_aligned_chains.cif'

    try:
        doc = open(alignedfile, 'rb')
    except: 
        print(f"Could find {filehandle}.Exited")
        return Response(-1)

    response = HttpResponse(FileWrapper(doc), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}-{}_{}-{}.cif"'.format(
    struct1,strand1,struct2,strand2)

    return response


    return Response("Hey now,\
    {} {} {} {}".format(
    handle1,handle2,
    struct1+"_"+strand1,
    struct2+"_"+strand2
    ))











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

    response = HttpResponse(FileWrapper(doc), content_type='application/json')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)

    return response



@api_view(['GET', 'POST'])
def tunnel(request):
    params     = dict(request.GET)
    struct     = params['struct'][0].upper()
    # centerline, report
    filetype   = params['filetype'][0]

    cetrline_filehandle  =  os.path.join(STATIC_ROOT, struct,'TUNNEL', 'csv', 'centerline.csv')
    report_filehandle    =  os.path.join(STATIC_ROOT, struct,'TUNNEL', f"{struct}_TUNNEL_REPORT.json")


    if filetype== 'report':
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



