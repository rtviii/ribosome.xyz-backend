import os
import requests

os.environ['NO_PROXY'] = '127.0.0.1'


# api-endpoint
URL = "http://localhost:8000/neo4j/cypher/"
# defining a params dict for the parameters to be sent to the API
PARAMS = {'cypher':"match (n:Protein) return n limit 1"}
# sending get request and saving the response as response object
r = requests.get(url = URL, params = PARAMS)
# extracting data in json format
data = r

print(data)