from neo4j import GraphDatabase, Result
import os

def _neoget(CYPHER_STRING:str):

    print("PASSSWORD IS ", os.getenv("NEO4J_PASSWORD"))
    driver = GraphDatabase.driver(
    #     os.getenv( 'NEO4J_URI' ),
    #     auth=(
    #     os.getenv( 'NEO4J_USER' ),
    #     os.getenv( 'NEO4J_PASSWORD')
    # )
    'bolt://ribosome.xyz/7687',auth=('rt','rrr')
    )

    def parametrized_query(tx, **kwargs):
        result:Result = tx.run(CYPHER_STRING, **kwargs)
        return result.values()
    with driver.session() as session:
        session.close()
        return session.read_transaction(parametrized_query)

def get_nom_cmap(rcsb_id:str)->dict:

    CYPHER_STRING = """
    match (n:RibosomeStructure{{rcsb_id:"{}"}})-[]-(r) where r:rRNA or r:RibosomalProtein with r 
    match (r)-[]-(d) where d:RPClass or d:RNAClass return r.entity_poly_strand_id, d.class_id""".format(rcsb_id)
    cmap:  dict  = {};
    response = _neoget(CYPHER_STRING)

    for resp_chain in response: 	
        if ','in resp_chain[0]:
            double_key = resp_chain[0].split(',')
            for key in double_key:
                cmap[key] = resp_chain[1]
                continue
        else:
            cmap[resp_chain[0]] = resp_chain[1]

    return cmap