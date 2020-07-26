from django.shortcuts import render
from rest_framework import viewsets          
from rest_framework.decorators import api_view
from rest_framework.response import Response
from .serializers import     MoleculeSerializer
from .models import Molecule     
from neo4j import GraphDatabase




class MoleculeView(viewsets.ModelViewSet):
    serializer_class = MoleculeSerializer
    queryset = Molecule.objects.all()


uri = "neo4j://localhost:7687"
driver = GraphDatabase.driver(uri, auth=("rt", "rrr"))





@api_view(['GET', 'POST'])
def get_molecule(request):
    print(request.GET)
    params = dict(request.GET)
    identifier = params['identifier']
    print("IDENTIFIER IS " ,identifier)
    CYPHER_STRING="""match (S:Subchain {nomenclature:%s})-[]-(pf:PfamFamily)
                with S, S.nomenclature as noms, pf
                return *;"""%identifier

    if request.method == 'GET':
        def get_molecules(tx):
            molecules = []
            return list(tx.run(CYPHER_STRING))
            
        with driver.session(database='ribosome') as session:
            result = session.read_transaction(get_molecules)
        driver.close()
        return Response(result)


