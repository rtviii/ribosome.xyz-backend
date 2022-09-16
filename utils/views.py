from rest_framework.decorators import api_view
from rest_framework.response import Response
from neo4j import  Result, GraphDatabase
import os
#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯
uri        =  os.getenv( 'NEO4J_URI'                                      )
authglobal = (os.getenv( 'NEO4J_USER'      ),os.getenv( 'NEO4J_PASSWORD' ))
current_db =  os.getenv( 'NEO4J_CURRENTDB'                                )

from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema

#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯

def _neoget(CYPHER_STRING:str):
    driver = GraphDatabase.driver(uri, auth= authglobal )
    def parametrized_query(tx, **kwargs):
        result:Result = tx.run(CYPHER_STRING, **kwargs)
        return result.value()
    with driver.session() as session:
        return session.read_transaction(parametrized_query)



@swagger_auto_schema(methods=[ 'get' ], auto_schema=None)  
@api_view(['GET'])
def number_of_structures(request):
    CYPHER_STRING="""match (r:RibosomeStructure) return count(r)"""
    return Response(_neoget(CYPHER_STRING))



