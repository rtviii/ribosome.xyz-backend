import collections
import functools
import json
import os
import numpy as np
from pymol import cmd, selector

# # TODO: -objectify singular chains,
#         -cli flag to highlight uL4 and uL22 and **ANY** other proteins?
#         -animate clusterwide-translation out along the reverse center-vector?

datapath         = "./clusterdata/"
nomenclatureMaps = './clusterdata/subchainMaps/'


def seeclusters(pdbid, radius=False):
    pdbid = str.upper(pdbid)
    radius = str(radius)
    nomMap = loadjson(nomenclatureMaps + "{}.json".format(pdbid))
    targetspath = datapath + "{}/".format(pdbid)

    available = os.listdir(targetspath)
    radNameMap = {}

    def extractRadius(filename):
        clusterdatum = loadjson(targetspath + filename)
        print("inspecting {}".format(filename))
        radstr = str(round(clusterdatum['metadata']['radius'], 2))
        radNameMap[radstr] = filename

    [extractRadius(x) for x in available]

    def availableSay():
        print("-----------------------------------")
        print("Avalaibale radii: \n")
        [print("{}".format(rad)) for rad in sorted(radNameMap.keys())]
        print("-----------------------------------")

    if radius not in radNameMap.keys():
        availableSay()
        print('Radius not found!')
        return

    report = loadjson(targetspath + radNameMap[radius])
    rnanames = report['metadata']['rnas']
    report_chains = report['metadata']['allchains']
    pymol_chains = [chain for chain in cmd.get_chains(pdbid)]

    print("# Chains:{} \n{} ".format(len(report_chains), report_chains))
    print("..of those --{} RNAs: {}".format(len(rnanames), rnanames))

    if len(pymol_chains) != len(report_chains):
        print("Report doesn't match RCSB-provided chains data. Exiting!")
        raise AssertionError

    pymol_proteins = [
        *filter(lambda chain: chain not in rnanames,  pymol_chains)]

    def chaintonomenclature(chainid):
        if chainid not in nomMap.keys():
            print("{} not found in nomenclature map".format(chainid))
            raise AttributeError
        if len(nomMap[chainid]) != 1:
            print("No identifier for {}".format(chainid))
            return []
        return nomMap[chainid]

    # pymol_ban = [* map(chaintonomenclature, pymol_proteins)]
    # print(pymol_ban)

    cmd.color("white", "all")
    for chain in rnanames:
        cmd.color('gray40', 'chain {}'.format(chain))
        cmd.hide("everything", "chain {}".format(chain))
        print("Hid rna {}".format(chain))

    def flattenDeep(l):
        for el in l:
            if isinstance(el, collections.Iterable) and not isinstance(el, (str, bytes)):
                yield from flattenDeep(el)
            else:
                yield el

    def getClusterChains(cluster):
        idscluster = []
        for nomid in cluster:
            # track back the nomenclature to the chain id
            if nomid not in flattenDeep([* nomMap.values()]):
                print("Ambiguous nomenclature {}".format(nomid))
                idscluster.append(nomid)
                pass
            for kvpair in [* nomMap.items()]:
                if nomid in kvpair[1]:
                    if len(kvpair[1]) > 1:
                        print("Ambiguous chain {} = {}".format(
                            nomid, kvpair[0]))
                    # if len(kvpair[1]) == 0:
                    #     idscluster.append()
                    idscluster.append(kvpair[0])
        return idscluster

    def revertNomenclature(chainid):
        if chainid not in flattenDeep([* nomMap.values()]):
            return chainid
        else:
            for kvp in [* nomMap.items()]:
                if chainid in kvp[1]:
                    return kvp[0]

    palette = ["purpleblue", "pink", "lightblue",  "greencyan", "blue", "palecyan",
               "violet",  "yelloworange", "limegreen", "warmpink", "raspberry", "lightmagenta",  "red", "palegreen", "wheat", "limon", "smudge"]

    # a book-keeping list to be emptied as clustered chains are converted into objects
    clustered_chains = []
    # for allchain in [* map(revertNomenclature, report['metadata']['allchains'])]:
    for allchain in report['metadata']['allchains']:
        if allchain in report['metadata']['singular_chains']:
            create_chain_object(revertNomenclature( allchain ), 'singular', nomMap)
        elif allchain in report['metadata']['rnas']:
            create_chain_object(revertNomenclature( allchain ), 'rna', nomMap)
        else:
            print("neither rna nor singular. adding to clusters {}".format(allchain))
            clustered_chains.append(allchain)
    
    print(report['clusters'])
    for clustpos, cluster in enumerate(report['clusters']):
        print("Cluster {}".format(cluster))
        for chain in cluster:
            create_chain_object(revertNomenclature( chain ), palette[clustpos], nomMap)
            print("attempting to remove {}".format(chain))
            clustered_chains.remove(chain)

    print(clustered_chains)
    print("why no loop")
    assert(len(clustered_chains) == 0)



# CHAINID = name of the chain as in pdb
# RNASINGCLUST = 'rna' | 'singular' | 'color'  if belongs to cluster associate with that color
# NOMENCLATUREMAP = object with chainid-bannoms relationships for the molecule(generated by the node module)
def create_chain_object(chain, rnasingclust, nomenclatureMap):
    nomMap = nomenclatureMap
    if rnasingclust != 'rna':
        ambiguous = True if len(nomMap[chain]) != 1 else False
    if rnasingclust == 'rna':
        cmd.select("rnasel_{}".format(chain), "chain {}".format(chain))
        cmd.create("rna_{}".format(chain), "rnasel_{}".format(chain))
        # cmd.delete("rnasel_{}".format(chain))
    elif rnasingclust == 'singular':
        if ambiguous:
            cmd.select("chain_{}".format(chain), "chain {}".format(chain))
            cmd.create("{}".format(chain), "chain_{}".format(chain))
            # cmd.delete("chain_{}".format(chain))
        else:
            cmd.select("chain_{}_{}".format(
                chain, nomMap[chain][0]), "chain {}".format(chain))
            cmd.create("{}_{}".format(chain, nomMap[chain][0]), "chain_{}_{}".format(
                chain, nomMap[chain][0]))
            # cmd.delete("chain_{}_{}".format(chain, nomMap[chain][0]))
    else:
        if ambiguous:
            cmd.select("chain_{}".format(chain), "chain {}".format(chain))
            cmd.create("{}_{}".format(rnasingclust, chain),
                       "chain_{}".format(chain))
            cmd.color(rnasingclust, "chain {}".format(chain))
            # cmd.delete("chain_{}".format(chain))
        else:
            cmd.select("chain_{}".format(chain), "chain {}".format(chain))
            cmd.create("{}_{}_{}".format(
                chain, nomMap[chain][0], rnasingclust), "chain_{}".format(chain))
            cmd.color(rnasingclust, "chain {}".format(chain))
            # cmd.delete("chain_{}".format(chain))


def objectify_chain_arr(chainid_arr=list):
    for chain in chainid_arr:
        cmd.select("chain_{}_selection".format(
            chain), "chain {}".format(chain))
        cmd.create("object_{}".format(chain),
                   "chain_{}_selection".format(chain))
        cmd.delete("chain_{}_selection".format(chain))


def loadjson(filename=str):
    with open(filename, 'r') as infile:
        jsonfile = json.load(infile)
    return jsonfile


def save_to_root(name):
    cmd.set("pse_export_version", 1.74)
    cmd.save(name + ".pse")


def whisp_menu():
    print("_____________________________________________________________________________________________________")
    print("|              SIGNATURE               |              EFFECT                                        |")
    print("|seeclusters [pdbid],                  | Highlight clusters, creating objects for chains of interest|")
    print("|save_under [filename],                | Save workspace with a name provided to the local folder    |")
    print("|___________________________________________________________________________________________________|")


cmd.extend("seeclusters", seeclusters)
cmd.extend("save_under", save_to_root)
whisp_menu()
