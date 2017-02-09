#!/home/dn070017/anaconda3/bin/python

import re
import sys

from collections import defaultdict

probe_file = open('/home/dn070017/projects/transcript_quatification/mouse_data/filter_genes_v1.fasta', 'r')
probe_key = ''
former_length = 0
probe_length = dict()
for probe_line in probe_file:
    probe_line = probe_line.strip()
    match_header = re.match('>(\S+)\|\S+', probe_line)
    if match_header:
        if probe_key != '':
            if probe_key not in probe_length:
                probe_length[probe_key] = former_length
                former_length = 0
            else:
                print('duplicate probe id: {}'.format(probe_key), file=sys.stderr)
        probe_key = match_header.group(1)
    else:
        former_length += len(probe_line)

probe_file.close()
if probe_key not in probe_length:
    probe_length[probe_key] = former_length
else:
    print('duplicate probe id: {}'.format(probe_key), file=sys.stderr)

blastn_file = open('/home/dn070017/projects/transcript_quatification/mouse_data/blastn_mouse_cdna_probe.tsv', 'r')
# query_choice = defaultdict(list)
reference_choice = defaultdict(list)
# scores = defaultdict(list)
for blastn_line in blastn_file:
    blastn_line = blastn_line.strip()
    blastn_data = blastn_line.split('\t')
    query = blastn_data[0]
    reference = blastn_data[1]
    identity = float(blastn_data[2])
    length = int(blastn_data[3])
    evalue = float(blastn_data[10])
    bitscore = float(blastn_data[11])

    if identity < 100 or evalue > 1e-5:
        continue

    query_split = query.split('|')
    query = query_split[0]
    
    if length < probe_length[query]:
        continue

    reference_split = reference.split('|')
    reference = reference_split[0] 

    reference_choice[reference].append(query)
    # print('{}\t\t{}\t{}\t{}\t{}\t{}'.format(query, reference, identity, length, evalue, bitscore))

blastn_file.close()

microarray_file = open('/home/dn070017/projects/transcript_quatification/mouse_data/mouse_real_microarray.tsv', 'r')
array_value_h = dict()
array_value_l = dict()
for microarray_line in microarray_file:
    microarray_line = microarray_line.strip()
    if microarray_line[0:9] == 'ProbeName':
        continue
    microarray_data = microarray_line.split('\t')
    array_value_h[microarray_data[0]] = float(microarray_data[3])
    array_value_l[microarray_data[0]] = float(microarray_data[4])
    # print('{}\t{}\t{}'.format(microarray_data[0], microarray_data[3], microarray_data[4]))

microarray_file.close()

print('target_id\tsample_01\tsample_2')
for reference, probe_list in reference_choice.items():
    array_value_h_sum = 0
    array_value_l_sum = 0
    for probe in probe_list:
        array_value_h_sum += array_value_h[probe]
        array_value_l_sum += array_value_l[probe]
    print('{}\t{:.3f}\t{:.3f}'.format(reference, array_value_h_sum/len(probe_list), array_value_l_sum/len(probe_list)))

