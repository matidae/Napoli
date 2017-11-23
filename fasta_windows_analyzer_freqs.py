import sys
import re
from math import ceil
from Bio import SeqIO

data = SeqIO.parse(sys.argv[1], 'fasta')
window_size = int(sys.argv[2])
print "\t".join(["Chromosome", "window","start","end", "A", "T", "C", "G", "N", "ATCG", "GC", "CG", "TG", "CTG", "CAG", "CCG", "CNG", "CXG", "CCGG", "CCCGGG"])
for i in data:
    name = i.id
    seq_len = len(i.seq)
    current_position = 0
    n_windows = int(ceil(seq_len*1.0/window_size))
    for j in xrange(n_windows):
        seq = str(i.seq[current_position : current_position + window_size]).upper()
        count_A = seq.count('A')
        count_T = seq.count('T')
        count_C = seq.count('C')
        count_G = seq.count('G')
        count_N = seq.count('N')
        count_ATCG = count_A + count_T + count_C + count_G
        count_GC = seq.count('GC')
        count_CG = seq.count('CG')
        count_TG = seq.count('TG')
        count_CTG = seq.count('CTG')
        count_CAG = seq.count('CAG')
        count_CCG = seq.count('CCG')
        count_CGG = seq.count('CGG')
        count_CNG = seq.count('CNG')
        count_CXG = len(re.findall('C.G', seq))
        count_CCGG = seq.count('CCGG')
        count_CCCGGG = seq.count('CCCGGG')
        if j != n_windows-1:
            end_position = current_position + window_size
        else:
            end_position = seq_len
        window_len = end_position - current_position
        counts_list =[count_A, count_T, count_C, count_G,\
                count_N, count_ATCG, count_GC, count_CG, count_TG, count_CTG,\
                count_CAG, count_CCG, count_CGG, count_CNG, count_CXG, count_CCGG,\
                count_CCCGGG]
        count_freqs = [h*1.0/window_len for h in counts_list]
        print("  ".join(map(str, [name, j+1, current_position+1, end_position] + count_freqs)))
        current_position += window_size
