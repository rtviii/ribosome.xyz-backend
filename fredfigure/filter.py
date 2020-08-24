import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

filename='./rcsb_pdb_custom_report_20200814185656.csv'
df = pd.read_csv(filename, usecols=['PDB ID', 'Deposition Date', 'Title', 'DOI'])
# df['Deposition Date']=pd.to_datetime(df['Deposition Date']) 

dois =df.groupby('DOI').agg({'PDB ID':lambda x: ' '.join(x), 'Deposition Date': lambda x: ' '.join(x)})
# dois['PDB ID'] =df['PDB ID'].apply(lambda x: list(dict.fromkeys( x.split(' ') )))

dois['PDB ID']= dois['PDB ID'].apply(lambda x: list(dict.fromkeys( x.split(' ') )))
dois['Struct Count'] = dois['PDB ID'].apply(lambda x: len(x))
dois['Deposition Date']  = dois['Deposition Date'].apply(lambda x: list( dict.fromkeys( map( lambda y: y.split('-')[0], x.split(' ') ) ) ))
dois =dois.sort_values(by="Struct Count", ascending=False)
print(dois)


dois.to_csv('byStructRefcount.csv')





