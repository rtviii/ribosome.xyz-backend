import os, sys,json,re
from dotenv import load_dotenv
import argparse


ROOT='/home/rxz/dev/ribxz'
sys.path.append(ROOT)

load_dotenv(os.path.join(ROOT,'.env'))

from ciftools.TunnelScripts.TunnelLog import ( Log)
from pymol import cmd

log=Log(os.getenv('TUNNEL_LOG'))

# Pymol's cmd-select should be a no-trace decorator. So sick of typing these cunts out
def tview(pdbid:str):

    species  = str( int( log.get_struct(pdbid)['taxid'].values[0] ) )

    l22chain = log.get_struct( pdbid ).uL22.values[0]
    l4chain  = log.get_struct( pdbid ).uL4.values[0]

    l22_res  = int( log.get_struct( pdbid ).constrictionResidueL22_id.values[0] )
    l4_res   = int( log.get_struct( pdbid ).constrictionResidueL4_id.values[0] )

    # cmd.delete('all')

    pdbid           = pdbid.upper()

    TUNNELS         = os.getenv("TUNNELS")
    SCOOP_RADIUS    = os.getenv("SCOOP_RADIUS")
    TUNNEL_SCRIPT   = os.path.join(TUNNELS,species, pdbid, 'pymol', 'complex.py')

    inputstructpath = os.path.join(TUNNELS,species, pdbid, '{}_{}Ascoop.pdb'.format(pdbid,SCOOP_RADIUS))
    csvpaths        = os.path.join(TUNNELS,species,pdbid,'csv',)





    cmd.load(inputstructpath)
    cmd.color('gray', 'all')

    if not os.path.exists(TUNNEL_SCRIPT):

        print("TUNNEL SCRIPT NOT FOUND.")
        return

    cmd.run(TUNNEL_SCRIPT)

    cmd.hide('everything', 'Tunnels')
    cmd.show('mesh', 'Tunnels')

    # tunnels    = os.listdir(csvpaths)
    # tunnumbers = map(lambda x: re.findall(r'\d+',x)[0], tunnels)
    # [cmd.delete(f'Tunnel{tunN}') if tunN not in choices else None for tunN in tunnumbers]
    # paint_tunnel(pdbid)


    sele_ptc(9606)
    # cmd.select('PTC', 'resi 2055 or resi 2056 or resi 2451 or resi 2452 or resi 2507 or resi 2506')
    # cmd.create('PTC',"PTC")
    # cmd.color('blue', 'PTC')

    cmd.select('ConstrictionSite', 'c. {} and resi {} or c. {} and resi {}'
    .format(l22chain,l22_res,l4chain,l4_res))

    cmd.create('ConstrictionSite',"ConstrictionSite")
    cmd.show('everything','ConstrictionSite')
    cmd.color('magenta', 'ConstrictionSite')

    cmd.reset()

def twrite(pdbid, args):
    # Anticipating the missing tunnels, usually would write PDBID, [1 2 3] to choose tunnel.
    # In the case of faulty struct/absence of tunnels : write PDBID, 0

    # Write -1 for a chain present inside the tunnel

    pdbid         = pdbid.upper()
    args          = args.split(' ')

    TUNNELS      = os.getenv("TUNNELS")
    CHOICE_TXT    = os.path.join(TUNNELS,species, pdbid,  'tunnels-results.txt')

    # if len(args) == 1 and int( args[0] ) == 0:
    # then it failed thats what's up. just write it

    f = open(CHOICE_TXT, 'w')
    f.write(",".join(args))
    f.close()
    print(f"Wrote results to {CHOICE_TXT}")
    cmd.delete('all')
    print("Cleared all, pymol.")

def choose_tunnel_mole(pdbid, args):


    log.update_struct(pdbid,moletunnel=[args])
    log._write()

    print("Wrote to {}".format(pdbid))

def paint_tunnel(pdbid:str):

    cmd.set('cartoon_transparency', '0.8','all')
    pdbid    = pdbid.upper()
    filename = '/home/rxz/dev/ribxz/static/{}/{}_TUNNEL_REPORT.json'.format(pdbid,pdbid)

    with open(filename) as report:
        report = json.load(report)


    cmd.select('sele_rna_lining',"c. A and resi 1000000")
    cmd.create('rna_lining', "sele_rna_lining")


    def paint_protein_chain(strand:str, banName:str):

        colormap ={
            'uL22': "gray70",
            'uL4' : 'gray70',
            'uL3' : "gray70",
            'uL23': "gray70"
        }

        if banName in colormap.keys():
            color = colormap[banName]

        else:
            color='pink'

        cmd.select("_sele_chain_{}".format(strand),
                   "c. {}".format(strand))

        cmd.create("{}__strand-{}".format(banName,strand), "_sele_chain_{}".format(strand),)

        cmd.color('gray70', "{}__strand-{}".format(banName,strand))
        cmd.hide('everything',"{}__strand-{}".format(banName,strand))
        cmd.show('sticks',"{}__strand-{}".format(banName,strand))

        residues = report['proteins'][strand]

        for res in residues:
            if res['resname'].upper() in ['LYS','ARG']:
                color ="red"
            cmd.select('c. {} and resi {}'.format(strand, res['resid']))
            cmd.color(color,'sele')

            if res['resname'].upper() in ['GLU','ASP']:
                color ="cyan"
            cmd.select('c. {} and resi {}'.format(strand, res['resid']))
            cmd.color(color,'sele')


    def paint_rna(strand:str):

        def update_rna_lining(strand,resid):
            cmd.select('sele_rna_lining','rna_lining or c. {} and resi {}'.format(strand, resid))
            cmd.create('rna_lining', 'sele_rna_lining')

        colormap={
            'A':'tv_green',
            'C':'tv_green',
            'G':'tv_green',
            'U':'tv_green'
        }

        for nucleotide in report['adjacent_strands'][strand]:

            if nucleotide['resname'] in colormap.keys():
                # target = 'sele_{}_{}'.format(strand,nucleotide['resid'] )
                # cmd.select(target, "c. {} and resi {}".format(strand, nucleotide['resid']))
                # cmd.color(colormap[nucleotide['resname']],target)
                # cmd.set('cartoon_transparency', '0.4',target)
                # cmd.delete(target)

                update_rna_lining(strand,nucleotide['resid'])
            else:
                target = 'sele_{}_{}'.format(strand,nucleotide['resid'] )
                cmd.select(target, "c. {} and resi {}".format(strand, nucleotide['resid']))
                cmd.create("{}".format(nucleotide['resname']), target)
                cmd.color("white",target)
                cmd.delete(target)

                

    rnaitems  = [* filter(lambda item: item[1]['type']=='RNA',report['nomenclatureMap'].items()) ]
    rnasnames = [*map(lambda item: item[0],rnaitems)]

    protitems = [* filter(lambda item: item[1]['type']=='Protein',report['nomenclatureMap'].items()) ]
    protnames = [*map(lambda item: item[0],rnaitems)]
    for rna in rnasnames:
        paint_rna(rna)

    for chain in protnames:
        nomenclature = report['nomenclatureMap'][chain]

        if len( nomenclature ) > 0:
            banName = report['nomenclatureMap'][chain]['nomenclature'][0]
    
        else:
            banName= None

        paint_protein_chain(chain,banName)


    cmd.hide('everything','rna_lining')
    cmd.show('surface','rna_lining')
    cmd.color('green', 'rna_lining')
    cmd.set('transparency','0.8', 'rna_lining')
    cmd.delete('sele')
    cmd.deselect()
    cmd.reset()


def sele_ptc(spec:int):
    # cmd.color('gray','all')
    if spec in [83333,562]:
        cmd.select('PTC', 'resi 2055 or resi 2056 or resi 2451 or resi 2452 or resi 2507 or resi 2506')
    if spec in [9606]:
        cmd.select('PTC', 'resi 4452')
    cmd.create('PTC',"PTC")
    cmd.color('blue', 'PTC')
    cmd.reset()

cmd.extend("tview",  tview)
cmd.extend("twrite", choose_tunnel_mole)
cmd.extend("tunnel", paint_tunnel)