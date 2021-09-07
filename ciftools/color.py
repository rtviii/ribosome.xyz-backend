import re


def hl_subseq(sequence:str, subsequence:str):
	CRED = '\033[91m'
	CEND = '\033[0m'
	_ = [ ]
	for item in re.split(re.compile(f'({subsequence})'),sequence):
		if item == subsequence:
			_.append(CRED + item + CEND)
		else:
			_.append(item)
	return ''.join(_)






seq1 = 'AAAAAAAFGJTTTKLCASCCCASCASCASIR1RFAS'
seq2 = 'TTT'

print(hl_subseq(seq1,seq2))