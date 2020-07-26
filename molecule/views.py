from django.shortcuts import render
from rest_framework import viewsets          
from rest_framework.decorators import api_view
from rest_framework.response import Response
from .serializers import     MoleculeSerializer
from .models import Molecule     
from neo4j import GraphDatabase
import os 




class MoleculeView(viewsets.ModelViewSet):
    serializer_class = MoleculeSerializer
    queryset = Molecule.objects.all()


driver = GraphDatabase.driver(os.environ['NEO4J_URL'],
 auth=(os.environ["NEO4J_USER"], os.environ['NEO4J_PASSWORD' ]))



@api_view(['GET', 'POST'])
def test_endp(request):
    return Response(os.environ)


@api_view(['GET', 'POST'])
def get_struct(request):
    params = dict(request.GET)
    pdbid = str.upper(params['pdbid'][0]) 

    
    CYPHER_STRING="""match (P:PDBStructure {{ pdbid:"{}"}})-[]-(s:Subchain)
                return * limit 25;""".format(pdbid)
    print(CYPHER_STRING)

    if request.method == 'GET':
        def make_query(tx):
            molecules = []
            return list(tx.run(CYPHER_STRING))
            
        with driver.session(database='ribosome') as session:
            result = session.read_transaction(make_query)
        driver.close()
        return Response(result)


