from neo4j import GraphDatabase
from rest_framework.decorators import api_view
from rest_framework.response import Response
import os
import environ
from Bio import PDB
# Neo4j driver has to be created inside the __init__ of the object to persist with its life. Otherwise is consumed on the first call.

environ.Env.read_env()

uri        = os.environ['NEO4J_URI']
authglobal = (os.environ['NEO4J_USER'],os.environ['NEO4J_PASSWORD'])
current_db = os.environ['NEO4J_CURRENTDB']

# driver = GraphDatabase.driver(uri,auth=('neo4j','55288'))

@api_view(['GET'])
def anything(request):
    print("OS environs :", os.environ )
    return Response('This is anything you want this to be;')

@api_view(['GET'])
def test_endpoint(request):
    return Response("This is testing endpoint. What did you expect?")

@api_view(['GET'])
def get_struct(request):
    params = dict(request.GET)
    pdbid  = str.upper(params['pdbid'][0])
    print("got params {}".format(params))

    # CYPHER_STRING="""match (P:PDBStructure {{ pdbid:"{}"}})-[]-(s:Subchain) return * limit 25;""".format(pdbid)
    CYPHER_STRING=f"""match(rr:rRNA)-[:rRNA_of]-(n:RibosomeStructure{{`_PDBId`:"{pdbid}"}})
                        with collect(rr) as rRNAs, n
                        match (rp:RibosomalProtein)-[:RibosomalProtein_of]-(n)
                        with collect(rp) as RPs, rRNAs, n
                        return {{ribosomalProteins: RPs,rRNAs: rRNAs,RibosomeStructure: n}};"""

    driver = GraphDatabase.driver(uri,auth=authglobal)
    def make_query(tx):
        molecules = []
        return list(tx.run(CYPHER_STRING))

    with driver.session(database=current_db) as session:
        result = session.read_transaction(make_query)
        driver.close()
    return Response(result)


@api_view(['GET'])
def get_homologs(request):
    params = dict(request.GET)
    ban = str(params['banName'][0])
    print('Got params {}'.format(params))
    CYPHER_STRING=f"""match (s:RibosomalProtein)-[]-(q:RibosomeStructure) where '{ban}' in s.nomenclature  
                    return {{protein:s, subchain_of: q.`_PDBId`}}"""
    driver = GraphDatabase.driver(uri,auth=authglobal)
    def make_query(tx):
        return tx.run(CYPHER_STRING)
    with driver.session(database=current_db) as session:
        result = session.read_transaction(make_query)
        driver.close()
    return Response(result)



@api_view(['GET','POST'])
def custom_cypher(request):
    driver = GraphDatabase.driver(uri,auth=authglobal)
    params = dict(request.GET)
    CYPHER_STRING  = params['cypher']

    def query(tx):
        results=[]
        return list(tx.run(CYPHER_STRING))
    with driver.session(database=current_db) as session:
        result = session.read_transaction(query)
    driver.close()
    return Response(result)

# @api_view(['GET'])
# def get_protein