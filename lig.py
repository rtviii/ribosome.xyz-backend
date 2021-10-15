import pandas as pd
import re

df = pd.read_csv('./unique_ligandlike.csv')


reg = re.compile(r"/(\w*(?<!(cha|pro|dom|str))in\b)|(\b\w*zyme\b)|(factor)/gi;")




def category(i:str)->str:

	if len(re.findall(r"(\w*(?<!(cha|pro|dom|str|pla))in\b|(\b\w*zyme\b))", i.lower()))> 0:
		return "Antibiotics"

	if len(re.findall(r"(factor)", i.lower()))> 0:
		return "Factors"

	if "mrna" in i.lower() or "messenger" in i.lower():
		return "mRNA"

	if "trna" in i.lower() or "transfer" in i.lower():
		return "tRNA"

	return "Mixed Ligands"


for idx,lig in enumerate(df.iloc[:,0]):
	df.iloc[idx,1] = category(lig)

df.to_csv('./unique_ligandlike.csv')