import re
import sys

from collections import defaultdict

if len(sys.argv) != 3:
    print('{} [trinity_nr.fasta] [answer_label.tsv]'.format(sys.argv[0]))
    sys.exit(0)

def cigar_to_pos(cigar_string):
    pos_list = list()
    curser = 0
    while True:
        cigar_match = re.match('(\d+)(\w)', cigar_string)
        if cigar_match:
            if cigar_match.group(2) == 'M':
                pos_list.append([curser + 1, curser + int(cigar_match.group(1))])
            else:
                curser += int(cigar_match.group(1))
        
            if cigar_match.end(0) == len(cigar_string):
                break
            cigar_string = cigar_string[cigar_match.end(0):]
    return pos_list

def check_overlap(a_start, a_end, b_start, b_end):
    a_length = a_end - a_start + 1
    b_length = b_end - b_start + 1
    average_length = (a_length + b_length) / 2
    if a_start >= b_start and a_end <= b_end:
        length = a_length
    if a_start >= b_start and a_end >= b_end:
        length = b_end - a_start + 1
    if a_start <= b_start and a_end >= b_end:
        length = b_length
    if a_start <= b_start and a_end <= b_end:
        length = a_end - b_start + 1
    if length / average_length < 0.1:
        return('separate')
    else:
        return('overlap')

def check_align_type(a_cigar, b_cigar):
    if a_cigar == b_cigar:
        return('overlap')
    
    a_pos_list = cigar_to_pos(a_cigar)
    b_pos_list = cigar_to_pos(b_cigar)

    for a_start, a_end in a_pos_list:
        for b_start, b_end in b_pos_list:
            overlap_type = check_overlap(a_start, a_end, b_start, b_end)
            if overlap_type == 'overlap':
                #print(a_pos_list, b_pos_list, overlap_type, (0, 0), sep='\t')
                return ('overlap', 0, 0)
    if a_end > b_start and a_start < b_start:
        #print(a_pos_list, b_pos_list, overlap_type, (int((a_end + b_start) / 2), -1), 'a', a_end, b_start, sep='\t')
        return ('separate', int((a_end + b_start) / 2), -1)
    if b_end > a_start and b_start < a_start:
        #print(a_pos_list, b_pos_list, overlap_type, (int((b_end + a_start) / 2), -1), 'b', b_end, a_start, sep='\t')
        return ('separate', int((b_end + a_start) / 2), -1)
    if a_end <= b_start and a_start < b_start:
        #print(a_pos_list, b_pos_list, overlap_type, (a_end, b_start), 'c', sep='\t')
        return ('separate', a_end, b_start)
    if b_end <= a_start and b_start < a_start:
        #print(a_pos_list, b_pos_list, overlap_type, (b_end, a_start), 'd', sep='\t')
        return ('separate', b_end, a_start)
    if a_start == b_start and a_end < b_end:
        #print(a_pos_list, b_pos_list, overlap_type, (a_end, -1), 'e', sep='\t')
        return ('separate', a_end, -1)
    if a_start == b_start and b_end < a_end:
        #print(a_pos_list, b_pos_list, overlap_type, (b_end, -1), 'e', sep='\t')
        return ('separate', b_end, -1)
    if a_end == b_end and a_start < b_start:
        #print(a_pos_list, b_pos_list, overlap_type, (b_start, -1), 'e', sep='\t')
        return ('separate', b_start, -1)
    if a_end == b_end and b_start < a_start:
        #print(a_pos_list, b_pos_list, overlap_type, (a_start, -1), 'e', sep='\t')
        return ('separate', a_start, -1)

def merge_junction(junctions):
    merged = list()
    for start, end in sorted(junctions):
        if start not in merged and start > 0:
            merged.append(start)
        if end not in merged and end > 0:
            merged.append(end)
    return(sorted(merged)) 

sequence_file = open(sys.argv[1], 'r')
seq_dict = dict()
header = ''
tmp_seq = ''
for sequence_line in sequence_file:
    sequence_line = sequence_line.strip()
    match = re.match('>(\S+)', sequence_line)
    if match:
        if header != '':
            seq_dict[header] = tmp_seq
        header = match.group(1)     
        tmp_seq = ''
    else:
        tmp_seq += sequence_line
sequence_file.close()

answer_label_file = open(sys.argv[2], 'r')
seq_length = dict()
target_cigar_list = defaultdict(list)
for answer_label_line in answer_label_file:
    answer_label_line = answer_label_line.strip()
    answer_label_data = answer_label_line.split('\t')
    if float(answer_label_data[7]) < 0.25 and float(answer_label_data[8]) < 0.25:
        continue
    ref = answer_label_data[0]
    target = answer_label_data[1]
    orientation = answer_label_data[2]
    seq_length[ref] = answer_label_data[3]
    seq_length[target] = answer_label_data[4]
    target_cigar_list[target].append((ref, answer_label_data[5], answer_label_data[6]))
answer_label_file.close()

align_type_dict = defaultdict(dict)
target_junctions = defaultdict(list)
for target, cigar_list in target_cigar_list.items():
    align_type_dict[target]['unique'] = 0
    align_type_dict[target]['overlap'] = 0
    align_type_dict[target]['separate'] = 0

    list_length = len(cigar_list)
    if list_length == 1:
        align_type_dict[target]['unique'] = 1
        continue

    for i in range(0, list_length):
        a_index = (target, cigar_list[i][0])
        for j in range(i + 1, list_length):
            b_index = (target, cigar_list[j][0])
            align_type = check_align_type(cigar_list[i][2], cigar_list[j][2])
            #align_type_dict[target][align_type[0]] += 1
            if align_type[0] == 'separate':
                tmp_junction = [align_type[1], align_type[2]]
                target_junctions[target].append(sorted(tmp_junction))
                #print(target, align_type[1], align_type[2])
      
    #print(target, align_type_dict[target], sep='\t')

merged_target_junctions = defaultdict(list)
for target, junctions in target_junctions.items():
    if len(junctions) != 1:
        continue
    merged = merge_junction(junctions)
    merged_target_junctions[target] = merged 

for header, sequence in seq_dict.items():
    if header not in merged_target_junctions:
        print('>' + header)
        print('\n'.join(sequence[i:i+80] for i in range(0, len(sequence), 80)))
    else:
        id = 1
        it = 0
        for junction in merged_target_junctions[header]:
            subsequence = sequence[it:junction]
            if len(subsequence) < 50:
                continue
            print('>' + header + '_split_' + str(id))
            print('\n'.join(subsequence[i:i+80] for i in range(0, len(subsequence), 80)))
            id += 1
            it = junction
        if it != len(sequence):
            subsequence = sequence[it:]
            if len(subsequence) < 50:
                continue
            print('>' + header + '_split_' + str(id))
            print('\n'.join(subsequence[i:i+80] for i in range(0, len(subsequence), 80)))

#print('target_id\tunique\tseparate\toverlap')
#for target, align_type in align_type_dict.items():
#    print(target, align_type['unique'], align_type['separate'], align_type['overlap'], sep='\t')
