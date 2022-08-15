from itertools import chain
from Bio import pairwise2
import gemmi
import pathlib
import argparse

# As per PDB 3J7Z ( https://www.rcsb.org/structure/3j7z )
ECOLI23SRRNA ="GGUUAAGCGACUAAGCGUACACGGUGGAUGCCCUGGCAGUCAGAGGCGAUGAAGGACGUGCUAAUCUGCGAUAAGCGUCGGUAAGGUGAUAUGAACCGUUAUAACCGGCGAUUUCCGAAUGGGGAAACCCAGUGUGUUUCGACACACUAUCAUUAACUGAAUCCAUAGGUUAAUGAGGCGAACCGGGGGAACUGAAACAUCUAAGUACCCCGAGGAAAAGAAAUCAACCGAGAUUCCCCCAGUAGCGGCGAGCGAACGGGGAGCAGCCCAGAGCCUGAAUCAGUGUGUGUGUUAGUGGAAGCGUCUGGAAAGGCGCGCGAUACAGGGUGACAGCCCCGUACACAAAAAUGCACAUGCUGUGAGCUCGAUGAGUAGGGCGGGACACGUGGUAUCCUGUCUGAAUAUGGGGGGACCAUCCUCCAAGGCUAAAUACUCCUGACUGACCGAUAGUGAACCAGUACCGUGAGGGAAAGGCGAAAAGAACCCCGGCGAGGGGAGUGAAAAAGAACCUGAAACCGUGUACGUACAAGCAGUGGGAGCACGCUUAGGCGUGUGACUGCGUACCUUUUGUAUAAUGGGUCAGCGACUUAUAUUCUGUAGCAAGGUUAACCGAAUAGGGGAGCCGAAGGGAAACCGAGUCUUAACUGGGCGUUAAGUUGCAGGGUAUAGACCCGAAACCCGGUGAUCUAGCCAUGGGCAGGUUGAAGGUUGGGUAACACUAACUGGAGGACCGAACCGACUAAUGUUGAAAAAUUAGCGGAUGACUUGUGGCUGGGGGUGAAAGGCCAAUCAAACCGGGAGAUAGCUGGUUCUCCCCGAAAGCUAUUUAGGUAGCGCCUCGUGAAUUCAUCUCCGGGGGUAGAGCACUGUUUCGGCAAGGGGGUCAUCCCGACUUACCAACCCGAUGCAAACUGCGAAUACCGGAGAAUGUUAUCACGGGAGACACACGGCGGGUGCUAACGUCCGUCGUGAAGAGGGAAACAACCCAGACCGCCAGCUAAGGUCCCAAAGUCAUGGUUAAGUGGGAAACGAUGUGGGAAGGCCCAGACAGCCAGGAUGUUGGCUUAGAAGCAGCCAUCAUUUAAAGAAAGCGUAAUAGCUCACUGGUCGAGUCGGCCUGCGCGGAAGAUGUAACGGGGCUAAACCAUGCACCGAAGCUGCGGCAGCGACGCUUAUGCGUUGUUGGGUAGGGGAGCGUUCUGUAAGCCUGCGAAGGUGUGCUGUGAGGCAUGCUGGAGGUAUCAGAAGUGCGAAUGCUGACAUAAGUAACGAUAAAGCGGGUGAAAAGCCCGCUCGCCGGAAGACCAAGGGUUCCUGUCCAACGUUAAUCGGGGCAGGGUGAGUCGACCCCUAAGGCGAGGCCGAAAGGCGUAGUCGAUGGGAAACAGGUUAAUAUUCCUGUACUUGGUGUUACUGCGAAGGGGGGACGGAGAAGGCUAUGUUGGCCGGGCGACGGUUGUCCCGGUUUAAGCGUGUAGGCUGGUUUUCCAGGCAAAUCCGGAAAAUCAAGGCUGAGGCGUGAUGACGAGGCACUACGGUGCUGAAGCAACAAAUGCCCUGCUUCCAGGAAAAGCCUCUAAGCAUCAGGUAACAUCAAAUCGUACCCCAAACCGACACAGGUGGUCAGGUAGAGAAUACCAAGGCGCUUGAGAGAACUCGGGUGAAGGAACUAGGCAAAAUGGUGCCGUAACUUCGGGAGAAGGCACGCUGAUAUGUAGGUGAGGUCCCUCGCGGAUGGAGCUGAAAUCAGUCGAAGAUACCAGCUGGCUGCAACUGUUUAUUAAAAACACAGCACUGUGCAAACACGAAAGUGGACGUAUACGGUGUGACGCCUGCCCGGUGCCGGAAGGUUAAUUGAUGGGGUUAGCGCAAGCGAAGCUCUUGAUCGAAGCCCCGGUAAACGGCGGCCGUAACUAUAACGGUCCUAAGGUAGCGAAAUUCCUUGUCGGGUAAGUUCCGACCUGCACGAAUGGCGUAAUGAUGGCCAGGCUGUCUCCACCCGAGACUCAGUGAAAUUGAACUCGCUGUGAAGAUGCAGUGUACCCGCGGCAAGACGGAAAGACCCCGUGAACCUUUACUAUAGCUUGACACUGAACAUUGAGCCUUGAUGUGUAGGAUAGGUGGGAGGCUUUGAAGUGUGGACGCCAGUCUGCAUGGAGCCGACCUUGAAAUACCACCCUUUAAUGUUUGAUGUUCUAACGUUGACCCGUAAUCCGGGUUGCGGACAGUGUCUGGUGGGUAGUUUGACUGGGGCGGUCUCCUCCUAAAGAGUAACGGAGGAGCACGAAGGUUGGCUAAUCCUGGUCGGACAUCAGGAGGUUAGUGCAAUGGCAUAAGCCAGCUUGACUGCGAGCGUGACGGCGCGAGCAGGUGCGAAAGCAGGUCAUAGUGAUCCGGUGGUUCUGAAUGGAAGGGCCAUCGCUCAACGGAUAAAAGGUACUCCGGGGAUAACAGGCUGAUACCGCCCAAGAGUUCAUAUCGACGGCGGUGUUUGGCACCUCGAUGUCGGCUCAUCACAUCCUGGGGCUGAAGUAGGUCCCAAGGGUAUGGCUGUUCGCCAUUUAAAGUGGUACGCGAGCUGGGUUUAGAACGUCGUGAGACAGUUCGGUCCCUAUCUGCCGUGGGCGCUGGAGAACUGAGGGGGGCUGCUCCUAGUACGAGAGGACCGGAGUGGACGCAUCACUGGUGUUCGGGUUGUCAUGCCAAUGGCACUGCCCGGUAGCUAAAUGCGGAAGAGAUAAGUGCUGAAAGCAUCUAAGCACGAAACUUGCCCCGAGAUGAGUUCUCCCUGACCCUUUAAGGGUCCUGAAGGAACGUUGAAGACGACGACGUUGAUAGGCCGGGUGUGUAAGCGCAGCGAUGCGUUGAGCUAACCGGUACUAAUGAACCGUGAGGCUUAACCU"
parser = argparse.ArgumentParser(description= 'CLI for locating PTC residues of 23SrRNA in a given prokaryotic PDB file')
parser.add_argument ("-t", "--targets", type= str, required=True)
parser.add_argument ("--display_all", type= bool)
args = parser .parse_args()
argdict    = vars(parser.parse_args())

if "targets" in argdict.keys():
    d["targets"] = [s.strip().upper() for s in d["targets"].split(",")]

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
            print("Found {} in {}".format(STRAND, chain_id))
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

    print("Got aligned IDs:",aln_ids)
    aln_ids = list(filter(lambda x: x != None, aln_ids ))

    for aln_resid in aln_ids:
        if tgt_aln[aln_resid] == '-':
            continue
        tgt_ids.append(backwards_match(tgt_aln,aln_resid))

    
    return [list(model[STRAND][ix][0].pos) for ix in tgt_ids]
    


for target in d["targets"]:
    target_ptc = process_target(target)
    print("[{}] PTC residues: {}".format(target,target_ptc))
    if args.display_all:
        print(target)
    print()

# cif = gemmi.cif.read_file('3J7Z.cif')
# block = cif.sole_block()

# --- alignment
# _            = pairwise2.align.globalxx(self.src,self.tgt, one_alignment_only=True)
# self.src_aln = _[0].seqA
# self.tgt_aln = _[0].seqB


# for ( chain_id,one_letter_code )  in zip(
#                                          block.find_loop('_entity_poly.pdbx_strand_id'),
#                                          block.find_loop('_entity_poly.pdbx_seq_one_letter_code')
#                                         ):
#     print(chain_id,one_letter_code)


# s1= "AAAAASDF"
# s2= "TAAAAASDF"
# s3= "TAAGGAASDF"

# s1 = SeqRecord(Seq(s1))
# s2 = SeqRecord(Seq(s2))
# s3 = SeqRecord(Seq(s3))

# maxlen = max(len(s1), len(s2), len(s3))
# for record in [s1,s2,s3]:
#     if len(record.seq) != maxlen:
#         sequence = str(record.seq).ljust(maxlen, '.')
#         record.seq = Seq(sequence)
# assert all(len(record.seq) == maxlen for record in [s1,s2,s3])

# ms= MultipleSeqAlignment([s1,s2,s3])

# for i in block:
#     print(block)


# with tempfile.NamedTemporaryFile() as outfile:
#     outfile.write(response.content)
#     outfile.flush()
#     filename = outfile.name
#     # cifstruct = FastMMCIFParser(QUIET=True).get_structure(RCSB_ID, filename)
#     # y = MMCIF2Dict.MMCIF2Dict(filename)
#     # doc   = gemmi.cif.read_file(filename)

#     cif_block = gemmi.cif.read(filename)[0]
#     structure = gemmi.make_structure_from_block(cif_block)
#     # for entity in structure.entities:
#     #     print(dir(entity))
#     #     print(entity.name)
#     #     print(entity.polymer_type)
#     #     print(entity.full_sequence)

#     for chain in structure:
#         print(chain)
