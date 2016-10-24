from collections import defaultdict
import re
import sys

if len(sys.argv) != 3:
    print('usage: {} [reference.fa] [blastn]'.format(sys.argv[0]))
    sys.exit(0)

name = ''
ref_length = dict()
ref_file = open(sys.argv[1], 'r')
for ref_line in ref_file:
    ref_line = ref_line.strip()
    id_match = re.match('>(\S+)', ref_line)
    if id_match:
        name = id_match.group(1)
        ref_length[name] = 0
    else:
        ref_length[name] += len(ref_line)
ref_file.close()

region_table = defaultdict(list)
blastn_file = open(sys.argv[2], 'r')
for blastn_line in blastn_file:
    blastn_line = blastn_line.strip()
    blastn_data = blastn_line.split('\t')
    name = blastn_data[0]
    start = int(blastn_data[6])
    end = int(blastn_data[7])

    if name not in region_table or [start, end] not in region_table[name]:
        region_table[name].append([start, end])

blastn_file.close()

print('target_id\tlength\talign_length\talign_coverage')
for name, length in ref_length.items():
    
    if name not in region_table:
        print(name, length, 0, 0.000, sep='\t')
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
    print(name, length, align_length, sep='\t', end='')
    print('\t{:.3f}'.format(align_length / length))
