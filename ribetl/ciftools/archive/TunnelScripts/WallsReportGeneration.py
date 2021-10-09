#!/usr/bin/env python3

import json
import os,sys
import pandas as pd
from dotenv import load_dotenv

def root_self(rootname:str='')->str:

    """Returns the rootpath for the project if it's unique in the current folder tree."""
    ROOT = os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

if __name__=="__main__":
    root_self('ribxz')


from ciftools.TunnelScripts.TunnelLog import (Log, TunnelWalls)
from ciftools.Structure import fetchStructure
from ciftools.Neoget import _neoget

# _562    = ["7K00", "5AFI", "3J9Y", "3J9Z", "3JA1", "3JCJ", "4UY8", "5GAD", "5GAE", "5GAG", "5GAH", "5JTE", "5JU8", "5L3P", "5LZA", "5LZD", "5LZE", "5MDV", "5MDW", "5MDZ", "5MGP", "5U9F", "5U9G", "5WDT", "5WE4", "5WE6", "5WF0", "5WFS", "6B4V", "6BOH", "6BOK", "6C4I", "6DNC", "6ENF", "6ENJ", "6ENU", "6GC0", "6GWT", "6GXM", "6GXN", "6GXO", "6HRM", "6I0Y", "6I7V", "6O9J", "6OGF", "6OGI", "6ORE", "6ORL", "6OSK", "6OSQ", "6OT3", "6OTR", "6OUO", "6OXA", "6OXI", "6PJ6", "6Q95", "6Q97", "6Q9A", "6S0K", "6U48", "6VU3", "6VWL", "6VWM", "6VWN", "6VYQ", "6VYR", "6VYS", "6WD0", "6WD1", "6WD2", "6WD3", "6WD4", "6WD5", "6WD7", "6WD8", "6WD9", "6WDA", "6WDD", "6WDE", "6WDK", "6WDM", "6WNT", "6WNV", "6WNW", "6X7F", "6X7K", "6XDQ", "6YSR", "6YSS", "6YST", "6YSU"]
# _83333  = ["3J7Z", "3JBU", "3JBV", "3JCD", "3JCE", "4U1U", "4U1V", "4U20", "4U24", "4U25", "4U26", "4U27", "4WF1", "4WOI", "4WWW", "4Y4O", "4YBB", "5CZP", "5DFE", "5FDU", "5FDV", "5GAG", "5H5U", "5IQR", "5IT8", "5J4D", "5J5B", "5J7L", "5J88", "5J8A", "5J91", "5JC9", "5JU8", "5KCR", "5KCS", "5KPS", "5KPW", "5KPX", "5L3P", "5MDY", "5NP6", "5NWY", "5O2R", "5U4I", "5U9F", "5U9G", "5UYK", "5UYL", "5UYM", "5UYP", "5UYQ", "5WFK", "6BU8", "6BY1", "6CFK", "6FKR", "6GBZ", "6GC8", "6OFX", "6OG7", "6SZS", "6WD6", "6WDF", "6WDG", "6WDJ", "6WDL", "6XZ7", "6XZA", "6XZB", "6Y69", "7BV8", "7JSS", "7JSW", "7JSZ", "7JT1", "7JT2", "7JT3"]
# _9606   = ["3J9M", "3J7Y", "3J92", "3JAG", "3JAH", "3JAI", "4UG0", "5A2Q", "5AJ0", "5LKS", "5LZT", "5LZU", "5LZV", "5LZW", "5LZX", "5LZY", "5LZZ", "5OOL", "5OOM", "5T2C", "6D90", "6D9J", "6EK0", "6G18", "6G5H", "6G5I", "6GAW", "6GAZ", "6GB2", "6HCF", "6I0Y", "6IP5", "6IP8", "6LQM", "6LSR", "6LSS", "6LU8", "6NU2", "6OLE", "6OLF", "6OLG", "6OLI", "6OLZ", "6OM0", "6OM7", "6P5N", "6QZP", "6R5Q", "6R6G", "6R6P", "6R7Q", "6RW4", "6RW5", "6T59", "6W6L", "6XA1", "6Y0G", "6Y2L", "6Y57", "6Z6L", "6Z6M", "6Z6N", "6ZVH", "7K5I"]
# _9986   = ["3J92", "3JAG", "3JAH", "3JAI", "3JAJ", "3JAN", "5FLX", "5LZS", "5LZT", "5LZU", "5LZV", "5LZW", "5LZX", "5LZY", "5LZZ", "6D90", "6D9J", "6GZ3", "6GZ4", "6GZ5", "6HCF", "6HCJ", "6MTB", "6MTC", "6MTD", "6MTE", "6P4G", "6P4H", "6P5I", "6P5J", "6P5K", "6P5N", "6R5Q", "6R6G", "6R6P", "6R7Q", "6SGC", "6T59", "6W2S", "6W2T"]
# _300852 = ["1VY4", "1VY5", "1VY6", "1VY7", "4P6F", "4P70", "4TUA", "4TUB", "4TUC", "4TUD", "4TUE", "4U1U", "4U1V", "4U20", "4U24", "4U25", "4U26", "4U27", "4V90", "4W29", "4W2E", "4W2F", "4W2G", "4W2H", "4W2I", "4W4G", "4WF1", "4WOI", "4WPO", "4WQ1", "4WQF", "4WQR", "4WQU", "4WQY", "4WR6", "4WRA", "4WRO", "4WSD", "4WSM", "4WT1", "4WT8", "4WU1", "4WZD", "4WZO", "4Y4O", "4Y4P", "4YPB", "4YZV", "4Z3S", "4Z8C", "4ZER", "4ZSN", "5CZP", "5DFE", "5DOX", "5DOY", "5E7K", "5E81", "5EL4", "5EL5", "5EL6", "5EL7", "5F8K", "5FDU", "5FDV", "5HAU", "5HCP", "5HCQ", "5HCR", "5HD1", "5IB8", "5IBB", "5IMQ", "5J30", "5J3C", "5J4B", "5J4C", "5J4D", "5J8B", "5MDY", "5NDJ", "5NDK", "5OT7", "5UQ7", "5UQ8", "5VP2", "5VPO", "5VPP", "5W4K", "5WIS", "5WIT", "5ZLU", "6BUW", "6BZ6", "6BZ7", "6BZ8", "6CAE", "6CFJ", "6CFK", "6CFL", "6CZR", "6FKR", "6GSJ", "6GSK", "6GSL", "6GZQ", "6ND5", "6ND6", "6NDK", "6NSH", "6NTA", "6NUO", "6NWY", "6O3M", "6O97", "6OF1", "6OF6", "6OJ2", "6OPE", "6ORD", "6OTR", "6OXA", "6OXI", "6Q95", "6QNQ", "6QNR", "6UCQ", "6UO1", "6XQD", "6XQE", "7JQL", "7JQM"]
# _559292 = ["3JCT", "4U3M", "4U3N", "4U3U", "4U4N", "4U4O", "4U4Q", "4U4R", "4U4U", "4U4Y", "4U4Z", "4U50", "4U51", "4U52", "4U53", "4U55", "4U56", "4U6F", "5APN", "5APO", "5DC3", "5FCI", "5FCJ", "5H4P", "5I4L", "5IT7", "5MC6", "5MEI", "5NDV", "5NDW", "5OBM", "5ON6", "5T62", "5TBW", "5Z3G", "6ELZ", "6EM1", "6EM3", "6EML", "6FAI", "6FYX", "6GQB", "6LQP", "6LQS", "6LQU", "6M62", "6N8J", "6N8K", "6N8L", "6N8M", "6N8N", "6N8O", "6RBD", "6RBE", "6Y7C", "6Z6J", "6Z6K", "6ZQB", "6ZQC", "6ZQD", "6ZQG", "7BT6", "7BTB"]


def InitWalls(pdbid:str)->TunnelWalls:
    """Initiate walls for a particular structure. This consumes a dataframe given that a choice of tunnel is present""" 

    def get_tunnels_dataframe(csvpath:str)->pd.DataFrame:
        tunnel_instance = pd.read_csv(csvpath)
        xyzr            = tunnel_instance[['Distance','FreeRadius', 'X','Y','Z']]
        return xyzr

    TUNNELS = os.getenv("TUNNELS")
    log     = Log(os.getenv('TUNNEL_LOG'))

    struct  = log.get_struct(pdbid)
    taxid   = str(int(struct['taxid'].values[0]))
    choice  = int(struct['moletunnel'].values[0] )

    tunnelfile = "tunnel_{}.csv".format(choice)

    if choice < 1:
        """See docs for disambiguation. Either mole hasn't succeeded or there is something blockign the tunnel."""
        print("Choice under 1.") 
        return

    tunnelcsv  = get_tunnels_dataframe(os.path.join(TUNNELS,taxid,pdbid,'csv',tunnelfile ))
    tw         = TunnelWalls(pdbid, fetchStructure(pdbid),tunnelcsv)
    tw.consumeMoleDataframe(10)
    return tw

def add_nomenclature_map_to_report(pdbid:str,path_to_report:str):

    with open(path_to_report, mode='rb') as reportfile:
        report = json.load(reportfile)

    chains =[]

    for chain in report['adjacent_strands'].keys():
        chains.append(chain)
    chains = ",".join(map(lambda x: "\"{}\"".format(x), chains ))

    nresponse = _neoget("""
    match (r:RibosomeStructure{{rcsb_id:"{}"}})-[]-(n) where n.entity_poly_strand_id in [{}]
    return {{strand:n.entity_poly_strand_id,type: n.entity_poly_polymer_type, nomenclature:n.nomenclature}};
    """.format(pdbid,chains))

    report['nomenclatureMap']={}

    for strand in nresponse:

        strand_profile = [*strand.values()][0]
        nomenclature   = strand_profile['nomenclature']
        strand         = strand_profile['strand']
        polytype       = strand_profile['type']

        report['nomenclatureMap'][strand] = {
            "type"        : polytype,
            "nomenclature": nomenclature
        }

    with open(path_to_report,'w') as reportfile:
        json.dump(report, reportfile)
    print("Added nomenclature map to {}".format(path_to_report))

# x = InitWalls(PDBID)