#!/home/dn070017/anaconda3/bin/python
import re
import sys

if len(sys.argv) != 4:
    print('usage: {} [sample_1.fasta] [sample_2.fasta] [reference.fasta]'.format(sys.argv[0]));
    sys.exit(0)

ref_file = open(sys.argv[3], 'r')
read_count_1 = dict()
read_count_2 = dict()
for ref_line in ref_file:
    ref_line = ref_line.strip()
    id_match = re.match('>(\S+)', ref_line)
    if id_match:
        read_count_1[id_match.group(1)] = 0
        read_count_2[id_match.group(1)] = 0
ref_file.close()

read_file = open(sys.argv[1], 'r')
for read_line in read_file:
    read_line = read_line.strip()
    id_match = re.match('>read\d+\/(\S+)', read_line)
    if id_match:
        read_count_1[id_match.group(1)] += 1
read_file.close()

read_file = open(sys.argv[2], 'r')
for read_line in read_file:
    read_line = read_line.strip()
    id_match = re.match('>read\d+\/(\S+)', read_line)
    if id_match:
        read_count_2[id_match.group(1)] += 1
read_file.close()

print('target_id\tsample_01\tsample_2')
for read_id, count in read_count_1.items():
    print('{}\t{}\t{}'.format(read_id, read_count_1[read_id], read_count_2[read_id]))
