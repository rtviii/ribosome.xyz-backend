#!/usr/bin/env python3
# This script is courtesy of the ribosome.xyz and its authors.
# This relise on the following packages to run
# - gemmi   : https://gemmi.readthedocs.io/en/latest/install.html
# - bipython: https://biopython.org/
# And additionally "requests" to download missing structures: https://pypi.org/project/requests/

# Distribute freely.

from Bio import pairwise2
import gemmi
import pathlib
import argparse

# As per PDB 3J7Z ( https://www.rcsb.org/structure/3j7z )
ECOLI23SRRNA ="GGUUAAGCGACUAAGCGUACACGGUGGAUGCCCUGGCAGUCAGAGGCGAUGAAGGACGUGCUAAUCUGCGAUAAGCGUCGGUAAGGUGAUAUGAACCGUUAUAACCGGCGAUUUCCGAAUGGGGAAACCCAGUGUGUUUCGACACACUAUCAUUAACUGAAUCCAUAGGUUAAUGAGGCGAACCGGGGGAACUGAAACAUCUAAGUACCCCGAGGAAAAGAAAUCAACCGAGAUUCCCCCAGUAGCGGCGAGCGAACGGGGAGCAGCCCAGAGCCUGAAUCAGUGUGUGUGUUAGUGGAAGCGUCUGGAAAGGCGCGCGAUACAGGGUGACAGCCCCGUACACAAAAAUGCACAUGCUGUGAGCUCGAUGAGUAGGGCGGGACACGUGGUAUCCUGUCUGAAUAUGGGGGGACCAUCCUCCAAGGCUAAAUACUCCUGACUGACCGAUAGUGAACCAGUACCGUGAGGGAAAGGCGAAAAGAACCCCGGCGAGGGGAGUGAAAAAGAACCUGAAACCGUGUACGUACAAGCAGUGGGAGCACGCUUAGGCGUGUGACUGCGUACCUUUUGUAUAAUGGGUCAGCGACUUAUAUUCUGUAGCAAGGUUAACCGAAUAGGGGAGCCGAAGGGAAACCGAGUCUUAACUGGGCGUUAAGUUGCAGGGUAUAGACCCGAAACCCGGUGAUCUAGCCAUGGGCAGGUUGAAGGUUGGGUAACACUAACUGGAGGACCGAACCGACUAAUGUUGAAAAAUUAGCGGAUGACUUGUGGCUGGGGGUGAAAGGCCAAUCAAACCGGGAGAUAGCUGGUUCUCCCCGAAAGCUAUUUAGGUAGCGCCUCGUGAAUUCAUCUCCGGGGGUAGAGCACUGUUUCGGCAAGGGGGUCAUCCCGACUUACCAACCCGAUGCAAACUGCGAAUACCGGAGAAUGUUAUCACGGGAGACACACGGCGGGUGCUAACGUCCGUCGUGAAGAGGGAAACAACCCAGACCGCCAGCUAAGGUCCCAAAGUCAUGGUUAAGUGGGAAACGAUGUGGGAAGGCCCAGACAGCCAGGAUGUUGGCUUAGAAGCAGCCAUCAUUUAAAGAAAGCGUAAUAGCUCACUGGUCGAGUCGGCCUGCGCGGAAGAUGUAACGGGGCUAAACCAUGCACCGAAGCUGCGGCAGCGACGCUUAUGCGUUGUUGGGUAGGGGAGCGUUCUGUAAGCCUGCGAAGGUGUGCUGUGAGGCAUGCUGGAGGUAUCAGAAGUGCGAAUGCUGACAUAAGUAACGAUAAAGCGGGUGAAAAGCCCGCUCGCCGGAAGACCAAGGGUUCCUGUCCAACGUUAAUCGGGGCAGGGUGAGUCGACCCCUAAGGCGAGGCCGAAAGGCGUAGUCGAUGGGAAACAGGUUAAUAUUCCUGUACUUGGUGUUACUGCGAAGGGGGGACGGAGAAGGCUAUGUUGGCCGGGCGACGGUUGUCCCGGUUUAAGCGUGUAGGCUGGUUUUCCAGGCAAAUCCGGAAAAUCAAGGCUGAGGCGUGAUGACGAGGCACUACGGUGCUGAAGCAACAAAUGCCCUGCUUCCAGGAAAAGCCUCUAAGCAUCAGGUAACAUCAAAUCGUACCCCAAACCGACACAGGUGGUCAGGUAGAGAAUACCAAGGCGCUUGAGAGAACUCGGGUGAAGGAACUAGGCAAAAUGGUGCCGUAACUUCGGGAGAAGGCACGCUGAUAUGUAGGUGAGGUCCCUCGCGGAUGGAGCUGAAAUCAGUCGAAGAUACCAGCUGGCUGCAACUGUUUAUUAAAAACACAGCACUGUGCAAACACGAAAGUGGACGUAUACGGUGUGACGCCUGCCCGGUGCCGGAAGGUUAAUUGAUGGGGUUAGCGCAAGCGAAGCUCUUGAUCGAAGCCCCGGUAAACGGCGGCCGUAACUAUAACGGUCCUAAGGUAGCGAAAUUCCUUGUCGGGUAAGUUCCGACCUGCACGAAUGGCGUAAUGAUGGCCAGGCUGUCUCCACCCGAGACUCAGUGAAAUUGAACUCGCUGUGAAGAUGCAGUGUACCCGCGGCAAGACGGAAAGACCCCGUGAACCUUUACUAUAGCUUGACACUGAACAUUGAGCCUUGAUGUGUAGGAUAGGUGGGAGGCUUUGAAGUGUGGACGCCAGUCUGCAUGGAGCCGACCUUGAAAUACCACCCUUUAAUGUUUGAUGUUCUAACGUUGACCCGUAAUCCGGGUUGCGGACAGUGUCUGGUGGGUAGUUUGACUGGGGCGGUCUCCUCCUAAAGAGUAACGGAGGAGCACGAAGGUUGGCUAAUCCUGGUCGGACAUCAGGAGGUUAGUGCAAUGGCAUAAGCCAGCUUGACUGCGAGCGUGACGGCGCGAGCAGGUGCGAAAGCAGGUCAUAGUGAUCCGGUGGUUCUGAAUGGAAGGGCCAUCGCUCAACGGAUAAAAGGUACUCCGGGGAUAACAGGCUGAUACCGCCCAAGAGUUCAUAUCGACGGCGGUGUUUGGCACCUCGAUGUCGGCUCAUCACAUCCUGGGGCUGAAGUAGGUCCCAAGGGUAUGGCUGUUCGCCAUUUAAAGUGGUACGCGAGCUGGGUUUAGAACGUCGUGAGACAGUUCGGUCCCUAUCUGCCGUGGGCGCUGGAGAACUGAGGGGGGCUGCUCCUAGUACGAGAGGACCGGAGUGGACGCAUCACUGGUGUUCGGGUUGUCAUGCCAAUGGCACUGCCCGGUAGCUAAAUGCGGAAGAGAUAAGUGCUGAAAGCAUCUAAGCACGAAACUUGCCCCGAGAUGAGUUCUCCCUGACCCUUUAAGGGUCCUGAAGGAACGUUGAAGACGACGACGUUGAUAGGCCGGGUGUGUAAGCGCAGCGAUGCGUUGAGCUAACCGGUACUAAUGAACCGUGAGGCUUAACCU"

parser = argparse.ArgumentParser(description= 'CLI for locating PTC residues of 23SrRNA in a given prokaryotic PDB file')
parser.add_argument ("-t", "--targets", type= str, required=True)
parser.add_argument ("--display_all", action='store_true')
args    = parser .parse_args()
argdict = vars(parser.parse_args())

if "targets" in argdict.keys():
    argdict["targets"] = [s.strip().upper() for s in argdict["targets"].split(",")]
    if len(argdict) > 50: 
        print("Please don't overload our servers. Paid out of pocket!:) \nInstead, get in touch for collaboration: rtkushner@gmail.com!")
        exit(1)

def backwards_match(alntgt:str, resid:int):
    """Returns the target-sequence index of a residue in the (aligned) target sequence"""
    if resid > len(alntgt):
        exit(IndexError(f"Passed residue with invalid index ({resid}) to back-match to target.Seqlen:{len(alntgt)}"))
    counter_proper = 0
    for i,char in enumerate(alntgt):
        if i == resid:
            return counter_proper
        if char =='-':
            continue
        else: 
            counter_proper  +=1

def forwards_match(alnsrc:str, resid:int):
    """Returns the index of a source-sequence residue in the aligned source sequence."""
    count_proper = 0
    for alignment_indx,char in enumerate( alnsrc ):
        if count_proper == resid:
            return alignment_indx
        if char =='-':
            continue
        else: 
            count_proper  +=1

def process_target(rcsb_id: str):
    default_path = f"{rcsb_id.upper()}.cif"
    if  not pathlib.Path(default_path).is_file():
        print(f"Could not locate file {default_path} in current directory. Downloading via {f'https://api.ribosome.xyz/static_files/download_structure?struct_id={rcsb_id}'}")
        import requests
        with open(default_path, 'wb') as outfile:
            outfile.write(requests.get(f'https://api.ribosome.xyz/static_files/download_structure?struct_id={rcsb_id}').content)

    target       = gemmi.cif.read_file(default_path)
    block        = target.sole_block()
    model        = gemmi.read_structure(default_path)[0]

    

    STRAND       = None
    SEQ          = None


    # Locate the chain of 23SrRNA class
    for (strand, nomclass) in zip(
        block.find_loop('_ribosome_nomenclature.entity_poly.pdbx_strand_id'),
        block.find_loop('_ribosome_nomenclature.polymer_class')
    ):
        if nomclass == '23SrRNA':
            STRAND = strand
            break

    # Now find sequence of this 23SrRNA
    for (chain_id, one_letter_code) in zip(
        block.find_loop('_entity_poly.pdbx_strand_id'),
        block.find_loop('_entity_poly.pdbx_seq_one_letter_code')
    ):
        if STRAND in chain_id.split(','):                      # X-RAY structures have 'dual' chains. Split on comma to check both.
            SEQ = str(one_letter_code).strip(";").strip("\n")

    if SEQ == None:
        print("Could not locate 23SrRNA sequence in {} CIF file".format(rcsb_id))

    alignment = pairwise2.align.globalxx(ECOLI23SRRNA,SEQ, one_alignment_only=True)
    src_aln      = alignment[0].seqA
    tgt_aln      = alignment[0].seqB
    
    ptc_ids = [2445,2446,2447,2448,2449,2450,2451,2452] # PTC residues in 23SrRNA
    aln_ids = []
    tgt_ids = []

    for src_resid in ptc_ids:
        aln_ids.append(forwards_match(src_aln,src_resid))
    aln_ids = list(filter(lambda x: x != None, aln_ids ))

    for aln_resid in aln_ids:
        if tgt_aln[aln_resid] == '-':
            continue
        tgt_ids.append(backwards_match(tgt_aln,aln_resid))

    
    return [list(model[STRAND][ix][0].pos) for ix in tgt_ids]
    

for target in argdict["targets"]:
    if not args.display_all:
        target_ptc = process_target(target)
        print("[\033[94m{}\033[0m] Approximate PTC position(1 of {} residues): \033[91m{}\033[0m".format(target,len( target_ptc ),target_ptc[0]))
    else:
        print("[\033[94m{}\033[0m] PTC atom positions: ".format(target))
        for residue in process_target(target):
            print(f"\t\033[91m{residue}\033[0m")
if not args.display_all:
    print("\n\nTo display more residues per target structure, use additional --display_all flag.")