from ast import operator
import json
from pprint import pprint
from unicodedata import name


with open('all_structs.json' , 'rb') as infile:
        data = json.load(infile)


def get_title(xobj):


        try:
                # return {xobj['data']['rcsb_id'] :xobj['data']['rcsb_primary_citation']['title']}
                return xobj['data']['rcsb_primary_citation']['title']
        except:
                return ""


# for i in data:
#         print(get_title(i))
titles =list(set([*map(lambda x: get_title(x), data)] ))







[print(y) for y in titles]










