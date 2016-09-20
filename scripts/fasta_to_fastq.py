#!/home/dn070017/anaconda3/bin/python
import os
import sys

if len(sys.argv) != 2:
    print('usage: {sys.argv[0]} [in.fasta]'.format(**locals()))
    sys.exit(0)

fasta_file_path = os.path.abspath(sys.argv[1])
fasta_file = open(fasta_file_path, 'r')
for fasta_line in fasta_file:
    fasta_line = fasta_line.strip()
    if fasta_line[0] == '>':
        fasta_line = '@' + fasta_line[1:]
    else:
        fasta_length = len(fasta_line)
        fasta_line = fasta_line + '\n+\n' + fasta_length * '?'
    print(fasta_line) 

fasta_file.close()
