import sys
import re

pat = re.compile("\s*\([0-9]*\) !")
subject = ""

with open(sys.argv[1]) as data:
    for line in data:
        if "Query:" in line:
            subject = line.split()[1]
        if re.match(pat, line):
            chromo = line.split()[5]
            if chromo != "cm" and chromo != "hmm":
                start = line.split()[6]
                end = line.split()[7]
                strand = line.split()[8]
                evalue = line.split()[2]
                score = line.split()[3]
                if strand == "+":
                    print("\t".join([chromo, "Infernal", "RNA", start, end, ".",\
                            strand, ".", "ID="+subject+";evalue="+evalue+";score="+score]))
                else:
                    print("\t".join([chromo, "Infernal", "RNA", end, start, ".",\
                            strand, ".", "ID="+subject+";evalue="+evalue+";score="+score]))
