#!/home/dn070017/anaconda3/bin/python
from collections import defaultdict
import re
import sys

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

if len(sys.argv) != 4:
    print('usage: {} [reference.fa] [blastn] [answer_label_details.tsv]'.format(sys.argv[0]))
    sys.exit(0)

details = set()
ambiguous_type = dict()
details_file = open(sys.argv[3], 'r')
for details_line in details_file:
    details_line = details_line.strip()
    details_data = details_line.split('\t')
    ambiguous_type[details_data[2]] = details_data[16]
    if details_data[16] != 'unique_ambiguous_ambiguous':
        continue
    details.add(details_data[2])

details_file.close()

name = ''
ref_seqs = dict()
ref_length = dict()
ambiguous_count = dict()
ref_file = open(sys.argv[1], 'r')
seq = ''
for ref_line in ref_file:
    ref_line = ref_line.strip()
    id_match = re.match('>(\S+)', ref_line)
    if id_match:
        if seq != '':
            ref_seqs[name] = seq
            seq = ''
        name = id_match.group(1)
        ref_length[name] = 0
        ambiguous_count[name] = 0
    else:
        seq += ref_line
        ref_length[name] += len(ref_line)
ref_seqs[name] = seq
ref_file.close()

ambiguous_type_count = dict()
blastn_file = open(sys.argv[2], 'r')
for blastn_line in blastn_file:
    blastn_line = blastn_line.strip()
    blastn_data = blastn_line.split('\t')

    reverse_complement = False

    name = blastn_data[0]
    t_name = blastn_data[1]
    identity = float(blastn_data[2])
    length = int(blastn_data[3])
    start = int(blastn_data[6])
    end = int(blastn_data[7])
    t_start = int(blastn_data[8])
    t_end = int(blastn_data[9])
    evalue = float(blastn_data[10])

    if name == t_name or length < 100:
        continue
    if evalue > 1e-5 or identity < 90:
        continue

    t_seq = ref_seqs[t_name]
    if t_start > t_end:
        reverse_complement = True
        original = 'ACGTacgt'
        convert = 'TGCAtgca'
        translator = str.maketrans(original, convert)
        t_start, t_end = t_end, t_start
        t_seq = t_seq[t_start:t_end]
        t_subseq = t_seq.translate(translator)[::-1]
    else:
        t_subseq = t_seq[t_start:t_end]
    
    ambiguous_count[name] += 1
    
    if name in details:
        print(name, t_name, ambiguous_type[t_name], reverse_complement, sep='\t')
        print('>' + ref_seqs[name][start:end], '>' + t_subseq, sep='\n')
        if ambiguous_type[t_name] not in ambiguous_type_count:
            ambiguous_type_count[ambiguous_type[t_name]] = 0
        else:
            ambiguous_type_count[ambiguous_type[t_name]] += 1

print(ambiguous_type_count, file=sys.stderr)
blastn_file.close()

"""
print('target_id\tlength\tambiguous_sites\talign_length\talign_coverage')
for name, length in ref_length.items():
    
    if name not in region_table:
        print(name, length, 0, 0, 0.000, sep='\t')
        continue

    regions = region_table[name]
    union_region = list()
    for start, end in sorted(regions):
        if union_region and union_region[-1][1] >= start - 1:
            union_region[-1][1] = max(union_region[-1][1], end)
        else:
            union_region.append([start, end])
        
    align_length = 0
    for start, end in union_region:
        align_length += (end - start + 1)
    
    print(name, length, ambiguous_count[name], align_length, sep='\t', end='')
    print('\t{:.3f}'.format(align_length / length))
"""        
