from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import MMCIF2Dict
from Bio.PDB.Structure import Structure
from Bio.PDB import FastMMCIFParser
from Bio.PDB import StructureBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Polypeptide import PPBuilder,Polypeptide
from rest_framework.decorators import api_view
from rest_framework.response import Response
import os

# @api_view(['GET'])
# def get_pdbsubchain(request):
#     pdbid = '3J9M'
#     filepath = f'./{pdbid}.cif'
#     cifparser = FastMMCIFParser(QUIET=True)
#     try:
#         with open(filepath) as infile:
#             structure:Structure = cifparser.get_structure(pdbid, infile)[0]
#     except:
#             print("Failed to open {}".format(filepath))

#     mychain:Chain = structure['D']

#     return Response('yep, you called it')


from rest_framework import generics
from django.http import HttpResponse
from wsgiref.util import FileWrapper

class FileDownloadListAPIView(generics.ListAPIView):

    def get(self, request, format=None, ):
        file_handle = './newchain.pdb'
        document = open(file_handle, 'rb')
        response = HttpResponse(FileWrapper(document), content_type='chemical/x-pdb')
        response['Content-Disposition'] = 'attachment; filename="newfile.pdb"'
        return response