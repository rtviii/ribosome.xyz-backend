import json
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib import rc, rcParams


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
    [year,month,day]  = datestr.split("-")
    if reso == None or int( reso[0] ) >4:
        print("fitlered out ", s)
        continue
    if int(year) < 2000:
        continue
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
rc('axes', linewidth=2)
# rc('font', weight='bold')
xraybars = ax.bar(d_years, d_xrays, width,                label='X-RAY DIFFRACTION' , fill = None   , edgecolor="black")
embars   = ax.bar(d_years, d_ems  , width, bottom=d_xrays,label='ELECTRON MICROSCOPY'  , color="black", edgecolor="black")
plt.xticks(d_years,fontsize=10)
# ax.set_ylabel('Scores')
# ax.set_title('Scores by group and gender')
ax.legend()

#  set x axis label to "Year"

ax.set_xlabel('Year'                                 , fontsize=16)
ax.set_ylabel('Number of entries'                    , fontsize=16)
title= "Number of entries in " + r"$\it{" + "ribosome.xyz" +"}$ by year"
ax.set_title (title, fontsize=18)
plt.legend(loc=2, prop={'size': 14})


# ax.set_xticklabels(['2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019'])
        

for ( bar, i_em, i_xray ) in zip(embars.patches, d_ems, d_xrays):
    value = bar.get_height()
    text = f'{i_em+i_xray}'
    text_x = bar.get_x() + bar.get_width() / 2
    text_y = bar.get_y() + value
    ax.text(text_x, text_y + 3, text, ha='center',color='black',size=12)

# plt.show()



fig.set_size_inches(18.5, 10.5)
fig.savefig('cumulative_entries.svg', format='svg', dpi=600)
fig.savefig('cumulative_entries.png', format='png', dpi=600)
fig.savefig('cumulative_entries.pdf', format='pdf', dpi=600)