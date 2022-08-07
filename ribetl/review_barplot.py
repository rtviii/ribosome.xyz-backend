import json
import numpy as np
import matplotlib.pyplot as plt
import datetime


# open the structs.json file and read the data 
with open('structs.json', 'r') as infile:
    data = json.load(infile)

dates   = []
heights = []
_       = 0

by_year  = {
    2000:{
        'xray': 0,
        'em'  : 0
    }
}


nmr  = 0
em   = 0
xray = 0


XRAY = 'X-RAY DIFFRACTION'
NMR  = 'SOLUTION NMR'
EM   = 'ELECTRON MICROSCOPY'


for s in data:
    _                += 1
    datestr           = str(s['data']['rcsb_accession_info']['deposit_date']).split("T")[0]
    method           = str(s['data']['exptl'][0]['method']); 
    if method == XRAY:
        xray +=1
    elif method == NMR:
        nmr +=1 
    elif method == EM:
        em+=1
        
        
    
    [year,month,day] = datestr.split("-")
    d                = datetime.datetime(int(year),int(month),int(day))
    dates.append(d)
    heights.append(_)
    



ax = plt.subplot(111)
ax.bar(dates, heights, width=10)
ax.xaxis_date()





    

# height = [3, 12, 5, 18, 45]
# bars   = ('A', 'B', 'C', 'D', 'E')
# x_pos  = np.arange(len(bars))
# plt.bar(x_pos, height, color=(0.1, 0.1, 0.1, 0.1),  edgecolor='blue')
# plt.xticks(x_pos, bars)
# plt.show()