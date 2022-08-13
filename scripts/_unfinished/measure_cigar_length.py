import sys

path_paf = sys.argv[1]

full_cigar_len_list = []
compressed_cigar_len_list = []

query_aligned_len_list = []
target_aligned_len_list = []

with open(path_paf) as f:
    for line in f:
        if 'cg:Z:' in line:
            query, query_len, query_start, query_end, query_strand, target, target_len, target_start, target_end = line.split('\t')[:9]

            full_cigar = line.split('cg:Z:')[1].split(' ')[0].strip()
            #print(full_cigar)

            # Compressed CIGAR
            cigar_tmp_list = []
            op_len = ''
            for x in full_cigar.replace('=', 'M').replace('X', 'M'):
                if x in '=XMID':
                    cigar_tmp_list.append((int(op_len), x))
                    op_len = ''
                else:
                    op_len += x
            compressed_cigar = ''
            prev_len, prev_op = cigar_tmp_list[0]
            for len_, op in cigar_tmp_list[1:]:
                if prev_op == op:
                    prev_len += len_
                else:
                    compressed_cigar += f'{prev_len}{prev_op}'
                    prev_len, prev_op = len_, op
            compressed_cigar += f'{prev_len}{prev_op}'

            query_aligned_len = int(query_end) - int(query_start)
            target_aligned_len = int(target_end) - int(target_start)
        else:
            full_cigar = ''
            compressed_cigar = ''
            query_aligned_len = 0
            target_aligned_len = 0

        full_cigar_len_list.append(len(full_cigar))
        compressed_cigar_len_list.append(len(compressed_cigar))

        query_aligned_len_list.append(query_aligned_len)
        target_aligned_len_list.append(target_aligned_len)

print('Average full CIGAR length: {:.3f}'.format(float(sum(full_cigar_len_list))/len(full_cigar_len_list)))
print('Average compressed CIGAR length: {:.3f}'.format(float(sum(compressed_cigar_len_list))/len(compressed_cigar_len_list)))

for i, (x, y, a, b) in enumerate(zip(full_cigar_len_list, compressed_cigar_len_list, query_aligned_len_list, target_aligned_len_list)):
    print(i+1, x, y, a, b)
