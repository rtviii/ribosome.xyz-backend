import pandas as pd
import numpy as np
from time import strptime
import re
pid='PDB ID'
mtd='Method'
yr='Year'
import matplotlib.pyplot as plt

# df = pd.read_csv('rcsb_pdb_custom_report_20200814185656.csv', usecols=['PDB ID', 'Deposition Date', 'Title', 'DOI', 'Assembly ID'])
# df['Deposition Date'] = df['Deposition Date'].apply(lambda datestring: str( strptime(datestring, "%Y-%m-%d").tm_year ))

# reduceOnId =df.groupby('PDB ID').agg({
#     'PDB ID'         : lambda x: ','.join(x),
#     'Deposition Date': lambda x: ','.join(x)
#     })


# reduceOnId['Deposition Date'] = reduceOnId['Deposition Date'].apply(lambda x: re.findall('\d{4}', x)[0])
# reduceOnId['PDB ID']          = reduceOnId['PDB ID'].apply(lambda x: re.findall('[A-Z|0-9]{4}',x)[0])
# reduceOnId.to_csv('reduceOnId.csv')



# PDB Requests found in expm.js


methods={
    "ELECTRON MICROSCOPY": 0,
    "X-RAY DIFFRACTION"  : 1
}

df = pd.read_csv('methods.csv', usecols=['PDB ID','Method','Year',])
df = df.sort_values(by='Year')

print(df)
tally = {}
for row in df.iterrows():
    experiment = row[1] 
    if experiment['Year'] not in tally.keys() and experiment['Year'] > 1999:
        tally[experiment['Year']] = [0,0]
    if experiment[mtd] not in methods.keys() or experiment['Year'] <2000 :
        print("Unrecognized: ", experiment[mtd])
        continue
    tally[experiment[yr]][methods[experiment[mtd]]] += 1;

tally = np.array([* map(lambda x: [x[0],*x[1]], tally.items()) ])
print(tally)
years = [*map(lambda x: int(x), tally[:,0])]
em    = [*map(lambda x: int(x), tally[:,1])]
xray  = [*map(lambda x: int(x), tally[:,2])]

all0 = 0;
em0 = 0;
for i in range(11):
    em0 += em[i]
    all0 += em[i] + xray[i]

all10 = 0;
em10 = 0;
for i in range(11,16):
    em10+= em[i]
    all10 += em[i] + xray[i]

all15 = 0;
em15 = 0;
for i in range(16,21):
    em15+= em[i]
    all15 += em[i] + xray[i]
print("\n\n")
print("2000 -- 2010 : {}%".format(round(100*  em0/all0 , 1)))
print("\t\t em : {} | total : {}".format(em0,all0))
print("2010 -- 2015 : {}%".format(round(100*  em10/all10, 1 )))
print("\t\t em : {} | total : {}".format(em10,all10))
print("2015 -- 2020 : {}%".format(round(100* em15/all15, 1)))
print("\t\t em : {} | total : {}".format(em15,all15))
print("\n\n")






# def autolabel(rects,annotations):
#     """Attach a text label above each bar in *rects*, displaying its height."""
#     for i,rect in enumerate(rects):
#         ax.annotate('{}'.format(annotations[i]),
#                     xy=(rect.get_x() + rect.get_width() / 2, annotations[i]),
#                     xytext=(0, 3),  # 3 points vertical offset
#                     fontsize=14,
#                     textcoords="offset points",
#                     ha='center', va='bottom')

# fig, ax = plt.subplots()
# w = 0.5
# x = np.arange(len(years))

# rects1 = ax.bar(x,height=xray,   width=w,  label='X-RAY DIFFRACTION', color='black',edgecolor='black')
# rects2 = ax.bar(x,    height=em,  width=w,bottom=xray, label='ELECTRON MICROSCOPY',color='white', edgecolor='black')
# autolabel( rects2,[ em[i] + xray[i] for i in range(len(xray)) ] )

# ax.set_ylabel('Number of Structures Deposited', fontsize=16)
# ax.set_xlabel('Year', fontsize=18)
# ax.set_xticks(x)
# ax.set_yticks(np.linspace(0,120,20,dtype=int))
# ax.set_yticklabels(np.linspace(0,120,20, dtype=int),fontsize=14)
# ax.set_xticklabels(years, fontsize=12)
# ax.legend(fontsize=18)
# ax.set_title("Structures Deposited by Method", fontsize=20)
# plt.show()





