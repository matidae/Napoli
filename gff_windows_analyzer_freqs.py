import sys
from math import ceil
from collections import OrderedDict

gff_file = open(sys.argv[1], 'r')
chromo_list = [(i.split()[1].rstrip(),int(i.split()[0])) for i in open(sys.argv[2])]
window_size = int(sys.argv[4])
genome = [[]] * len(chromo_list)
genome_win = [[]] * len(chromo_list)
repeat_dict = {i.split()[0]:i.split()[1] for i in open(sys.argv[3])} 
repeat_counter = {i.split()[0]:0 for i in open(sys.argv[3])} 

for i in xrange(len(chromo_list)):
    chromo_n_window = int(ceil(chromo_list[i][1] * 1.0 / window_size))
    genome[i] = [{val:0 for val in repeat_dict.values()} for _ in xrange(chromo_n_window)]
    start_win = 1
    for j in xrange(chromo_n_window):
        if j < chromo_n_window -1:
            genome_win[i].append([start_win, start_win + window_size -1])
            start_win += window_size
        else:
            genome_win[i].append([start_win, chromo_list[i][1]])

for i in gff_file:
    chromo_name = i.split()[0]
    chromo_number = int(chromo_name[-2:])
    start = int(i.split()[3])
    end = int(i.split()[4])
    window_number_start = (start/window_size) + 1
    window_number_end = (end/window_size) + 1
    desc = i.split()[8].split(";")[0].split("=")[1].split(":")[0]
    if desc in repeat_dict.keys():
        if window_number_start == window_number_end:
            genome[chromo_number-1][window_number_start-1][repeat_dict[desc]] += end - start + 1
        else:
            max_window = genome_win[chromo_number-1][window_number_start-1][1]
            genome[chromo_number-1][window_number_start-1][repeat_dict[desc]] += max_window - start + 1
            genome[chromo_number-1][window_number_end-1][repeat_dict[desc]] += end - max_window + 1

c = 0
#od_repeat = OrderedDict(sorted(repeat_dict.values()))
od_repeat = sorted(set([i[1] for i in list(repeat_dict.items())]))
print("\t".join(['Chromosome', 'window','start', 'end'] + list(od_repeat)))
for chromo in genome:
    w = 1
    for window in chromo:
        num = int(chromo_list[c][0][-2:])-1
        od_counts = OrderedDict(sorted(window.items()))
#        print od_counts
        od_freqs = ["{0:.4f}".format(val*1.0/(genome_win[num][w-1][1]-genome_win[num][w-1][0]+1)) for val in od_counts.values()]
        print("\t".join(map(str,[chromo_list[c][0], w, genome_win[num][w-1][0], genome_win[num][w-1][1]])) +"\t" + "\t".join(map(str,od_freqs))) 
        w += 1
    c += 1
    
