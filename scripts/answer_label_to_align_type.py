import re
import sys

from collections import defaultdict

if len(sys.argv) != 2:
    print('{} [answer_label.tsv]'.format(sys.argv[0]))
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
                return 'overlap' 
    return 'separate'

answer_label_file = open(sys.argv[1], 'r')
seq_length = dict()
target_cigar_list = defaultdict(list)
for answer_label_line in answer_label_file:
    answer_label_line = answer_label_line.strip()
    answer_label_data = answer_label_line.split('\t')
    ref = answer_label_data[0]
    target = answer_label_data[1]
    orientation = answer_label_data[2]
    seq_length[ref] = answer_label_data[3]
    seq_length[target] = answer_label_data[4]
    target_cigar_list[target].append((ref, answer_label_data[5], answer_label_data[6]))
answer_label_file.close()

align_type_dict = defaultdict(dict)
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
            align_type_dict[target][align_type] += 1
      
    #print(target, align_type_dict[target], sep='\t')

print('target_id\tunique\tseparate\toverlap')
for target, align_type in align_type_dict.items():
    print(target, align_type['unique'], align_type['separate'], align_type['overlap'], sep='\t')
