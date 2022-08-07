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
    '2000':{
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
    method            = str(s['data']['exptl'][0]['method']);
    reso = s['data']['rcsb_entry_info']['resolution_combined']
    if reso == None or int( reso[0] ) >4:
        continue
    print(reso)
    [year,month,day]  = datestr.split("-")
    d                 = datetime.datetime(int(year),int(month),int(day))
    print(reso)
    if year not in by_year:
        by_year[year] = {
            'xray': 0,
            'em'  : 0
        }

    if method == XRAY:
        by_year[year]['xray'] += 1
        xray +=1

    elif method == NMR:
        nmr +=1 
        
    elif method == EM:
        by_year[year]['em'] += 1
        em+=1
    
    dates.append(d)
    heights.append(_)
    
by_year = sorted(by_year.items())
# ------------------------------------------------------------

d_years = []
d_ems   = []
d_xrays = []

for e in by_year:
    _em   = e[1]['em']
    _xray = e[1]['xray']
    d_years.append(int(e[0]))
    d_ems.append(_em)
    d_xrays.append(_xray)

width         = 0.35       # the width of the bars: can also be len(x) sequence
fig, ax = plt.subplots()

xraybars = ax.bar(d_years, d_xrays, width,                label='xray' , fill = None   , edgecolor="black")
embars   = ax.bar(d_years, d_ems  , width, bottom=d_xrays,label='ems'  , color="black", edgecolor="black")
# ax.set_ylabel('Scores')
# ax.set_title('Scores by group and gender')
ax.legend()


        

for bar in embars.patches:
    value = bar.get_height()
    text = f'{value}'
    text_x = bar.get_x() + bar.get_width() / 2
    text_y = bar.get_y() + value
    ax.text(text_x, text_y + 3, text, ha='center',color='r',size=12)


plt.show()
