from Bio import SeqIO
from pprint import pprint


recordd = SeqIO.to_dict(SeqIO.parse('uL22_ALGN.fasta','fasta'))
print(recordd['5JC9'].seq)
print(recordd['5T2A'].seq)

print("**************")

recordd = SeqIO.to_dict(SeqIO.parse('uL4_ALGN.fasta','fasta'))
print(recordd['5JC9'].seq)
print(recordd['5T2A'].seq)