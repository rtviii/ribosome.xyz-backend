# from neo4j import GraphDatabase

# NEO4J_URI      = "bolt://ribosome.xyz:7687/"
# NEO4J_USER     = "rt"
# NEO4J_PASSWORD = "rrr"



# driver = GraphDatabase.driver( NEO4J_URI ,auth=(NEO4J_USER,NEO4J_PASSWORD))

# def get_classes(pdbid:str):
# 	def parametrized_query(tx, **kwargs):
# 		result = tx.run(
# 			"""match (r:RibosomeStructure{rcsb_id: "3J7Z"})-[]-(rp:RibosomalProtein)-[]-(n:NomenclatureClass) 
# 				return n.class_id, rp.entity_poly_strand_id
# 				""", **kwargs)
# 		return result.values()

# 	with driver.session() as session:
# 		session.close()

# 		return session.read_transaction(parametrized_query)


# def get_rna(pdbid:str):
# 	def parametrized_query(tx, **kwargs):
# 		result = tx.run(
# 				f"""match (n:rRNA)-[]-(rib:RibosomeStructure{{rcsb_id:'{pdbid}'}}) 
# 				return {{
# 					pdbx_description:n.rcsb_pdbx_description,
# 					strand_id       :n.entity_poly_strand_id}}"""				, **kwargs)
# 		return result.values()

# 	with driver.session() as session:
# 		session.close()
# 		return session.read_transaction(parametrized_query)
        