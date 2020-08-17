#!/usr/bin/env python
import json
import os
import sys
import argparse
import xml
from pymol import cmd
from makeconfig import make_input_config
import asyncio


sys.path.append(os.path.join(sys.path[0], 'mole2'))
sys.path.append(sys.path[0])
path_moleexe = os.path.join(sys.path[0], 'mole2', 'mole2.exe')

# path_moleconfig = os.path.join(sys.path[0],'input.xml')
path_moleconfig = os.path.join(sys.path[0], 'input-multiorigin.xml')
path_pdbstruct = os.path.join(os.path.dirname(
    sys.path[0]), 'static', 'pdb-structs', 'radiusobject.pdb')


def makeparxser():
    def ResidueSpecification(string: str):
        try:
            chainId, seqnum = map(str, string.split(','))
            return chainId, seqnum
        except:
            print("Malformed or absent origin Residue specification")
            raise SyntaxError

    parser = argparse.ArgumentParser(
        "See MOLE2.pdf documentation for detailed options")

    # Paths
    parser.add_argument('-output', '--output_path',dest='Output',
                        help='output folder for mole')
    parser.add_argument('-input', '--input_path', dest='Input',
                        help='input file to operate on')
    # Params ----
    # Cavity
    parser.add_argument('-pr', dest='ProbeRadius', nargs='?', required=False)
    parser.add_argument('-it', dest='InteriorThreshold',
                        nargs='?', required=False)
    parser.add_argument('-md', dest='MinDepth', nargs='?', required=False)
    # Tunnel
    parser.add_argument('-or', dest='OriginRadius', nargs='?', required=False)
    parser.add_argument('-scr', dest='SurfaceCoverRadius',
                        nargs='?', required=False)
    parser.add_argument('-br', dest='BottleneckRadius',
                        nargs='?', required=False)
    parser.add_argument('-mtl', dest='MinTunnelLength',
                        nargs='?', required=False)
    parser.add_argument('-mpl', dest='MinPoreLength',
                        nargs='?', required=False)
    parser.add_argument('-mts', dest='MaxTunnelSimilarity',
                        nargs='?', required=False)
    parser.add_argument('-bt', dest='BottleneckTolerance',
                        nargs='?', required=False)

    # Export
    parser.add_argument('--exports', dest='exports', choices=['t', 'tc', 'tcp'],
                        help="""
    Specify exports to enable thus
    [flag -> mole option]:
    t  -> tunnels only
    tc -> tunnesl and cavities
    tcp -> tunnels, cavities, poresauto, poresuser, poresmerged
      """)

    # Origins
    parser.add_argument('-o_points', dest='Points', nargs='*', action='append',
                        help='array of points to use as origin',
                        required=False)
    parser.add_argument('-o_residues', dest='Residues', nargs='*', action='append', type=ResidueSpecification,
                        help="""pass a comma-delimited string of Chain-name and Sequence number in that chain for each residue
    Ex. ...-o_residues A,130 A,145 A,160 ... results in <Residue Chain="A" SequenceNumber="130">\n
    <Residue Chain="A" SequenceNumber="145">\n 
    <Residue Chain="A" SequenceNumber="160">etc.
    """,
                        required=False)


    return parser


if __name__ == "__main__":
    moleparser = makeparxser()
    args       = moleparser.parse_args()
    #get rid of the unutilized arg options
    args       = filter((lambda kvpair: None not in kvpair), vars(args).items())
    args       = dict(args)
    # if 'help' in args:
    #     pass
        # optionsHelp()
    asyncio.run(make_input_config(args))
    os.system("mono ./mole/mole2/mole2.exe ./input-auto.xml")
