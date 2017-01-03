#!/home/dn070017/anaconda3/bin/python
from collections import defaultdict, OrderedDict
import math
import re
import sys

def usage():
    if len(sys.argv) != 4:
        print('usage: {} [reference.fa] [query.fa] [blastn]'.format(sys.argv[0]))
        sys.exit(0)
    return

def extract_sequence_length(file_path):
    name = ''
    fasta_length = dict()
    fasta_file = open(file_path, 'r')
    for fasta_line in fasta_file:
        fasta_line = fasta_line.strip()
        id_match = re.match('>(\S+)', fasta_line)
        if id_match:
            name = id_match.group(1)
            fasta_length[name] = 0
        else:
            fasta_length[name] += len(fasta_line)
    fasta_file.close()
    return fasta_length

def parsing_blast(file_path):
    best_evalue = dict()
    best_orientation = dict()
    q_region_table = defaultdict(list)
    r_region_table = defaultdict(list)
    blastx_file = open(file_path, 'r')
    for blastx_line in blastx_file:
        blastx_line = blastx_line.strip()
        blastx_data = blastx_line.split('\t')
        
        identity = float(blastx_data[2])
        evalue = float(blastx_data[10])
       
        if identity < 90 or evalue > 1e-5:
            continue
 
        q_name = blastx_data[0]
        r_name = blastx_data[1]

        q_start = int(blastx_data[6])
        q_end = int(blastx_data[7])

        r_start = int(blastx_data[8])
        r_end = int(blastx_data[9])

        orientation = 'forward'

        if r_start > r_end:
            r_start, r_end = r_end, r_start
            orientation = 'reverse'
        
        q_key = (q_name, r_name, orientation)
        r_key = (r_name, q_name, orientation)

        if q_name not in best_evalue or evalue < best_evalue[q_name]:
            best_orientation[q_name] = orientation
            best_evalue[q_name] = evalue
        
        if q_key not in q_region_table or [q_start, q_end] not in q_region_table[q_key]:
            q_region_table[q_key].append([q_start, q_end])
        if r_key not in r_region_table or [r_start, r_end] not in r_region_table[r_key]:
            r_region_table[r_key].append([r_start, r_end])

    blastx_file.close()
    return q_region_table, r_region_table, best_orientation

def union_region(region):
    union = list()
    for start, end in sorted(region):
        if union and union[-1][1] >= start - 1:
            union[-1][1] = max(union[-1][1], end)
        else:
            union.append([start, end])
    
    return union

def region_to_string_and_coverage(union, length):
    curser = 1
    align_num = 0
    align_string = ''
    for start, end in union:
        if start > curser:
            align_string += str(start - curser) + 'I'
        align_string += str(end - start + 1) + 'M'
        align_num += end - start + 1
        curser = end + 1
    if length >= curser:
        align_string += str(length - curser + 1) + 'I'
    
    return align_string, float(align_num/length) 


def combine_region(q_region_table, r_region_table, q_seq_length, r_seq_length, best_orientation):
    for key, q_regions in OrderedDict(sorted(q_region_table.items())).items():
        q_name = key[0]
        r_name = key[1]
        orientation = key[2]

        if best_orientation[q_name] != orientation:
            continue

        q_len = q_seq_length[q_name]
        r_len = r_seq_length[r_name]

        r_regions = r_region_table[(r_name, q_name, orientation)]

        q_union = union_region(q_regions)
        r_union = union_region(r_regions)

        q_string, q_coverage = region_to_string_and_coverage(q_union, q_len)
        r_string, r_coverage = region_to_string_and_coverage(r_union, r_len)
        
        print(q_name, r_name, orientation, q_len, r_len, q_string, r_string, sep='\t', end='\t')
        print('{:.3f}\t{:.3f}'.format(q_coverage, r_coverage))

    return

usage()
reference_length = extract_sequence_length(sys.argv[1])
query_length = extract_sequence_length(sys.argv[2])
q_region_table, r_region_table, best_orientation = parsing_blast(sys.argv[3])
combine_region(q_region_table, r_region_table, query_length, reference_length, best_orientation)

