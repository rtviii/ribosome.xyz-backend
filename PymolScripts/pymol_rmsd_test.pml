reinitialize

fetch 1F9J, async=0
fetch 1YX5, async=0

extract 1F9J_A, 1F9J and chain A
extract 1YX5_B, 1YX5 and chain B

test=cmd.super("1F9J_A","1YX5_B")

python
writefile=open("rmsd_file.txt","a")
writefile.write(' '.join('%s' % x for x in test))
writefile.write('\n')
writefile.close()
python end