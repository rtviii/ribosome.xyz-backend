from rest_framework.decorators import api_view
from neo4j import GraphDatabase
from rest_framework.response import Response
from mole.call_mono_mole import command
import os


driver = GraphDatabase.driver(os.environ['NEO4J_URL'],
 auth=(os.environ["NEO4J_USER"], os.environ['NEO4J_PASSWORD' ]))


@api_view(['GET'])
def test_endpoint(request):

    # call_mono_mole()
    os.system(command)

    return Response("This is testing endpoint. What did you expect?")

@api_view(['GET'])
def get_structs(request):
    params = dict(request.GET)
    pdbid  = str.upper(params['pdbid'][0])
    CYPHER_STRING="""match (P:PDBStructure {{ pdbid:"{}"}})-[]-(s:Subchain)
                return * limit 25;""".format(pdbid)
    def make_query(tx):
        molecules = []
        return list(tx.run(CYPHER_STRING))
    with driver.session(database='ribosome') as session:
        result = session.read_transaction(make_query)
    driver.close()
    return Response(result)

@api_view(['GET'])
def get_all_outgoing_struct(request):

    params = dict(request.GET)
    pdbid  = str.upper(params['pdbid'][0])

    CYPHER_STRING="""match (p:PDBStructure{{pdbid:"{}"}})-[issubchain:Is_Subchain_Of]-(chain:Subchain)-[belong:Belongs_to]-(pfam:PfamFamily)-[ifam:InterProFamily_PfamFamily_CrossRef]-(ipro:InterProFamily)-[igo:InterProFamily_GoClass_CrossRef]-(goclass:GoClass)
    with p,issubchain, chain, belong, pfam,ifam,ipro,igo,goclass
    return * limit 400;""".format(pdbid)
    def query(tx):
        results=[]
        return list(tx.run(CYPHER_STRING))
    with driver.session(database='ribosome') as session:
        result = session.read_transaction(query)
    driver.close()
    return Response(result)