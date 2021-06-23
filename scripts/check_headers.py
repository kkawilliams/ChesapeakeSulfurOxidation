import os
import sys
from Bio import SeqIO

args = sys.argv
files_ = [os.path.join(args[1], i) for i in os.listdir(args[1]) if i.endswith(".fastq.test")]
assert len(files_) == 3
assert sum([os.path.exists(i) for i in files_]) == 3
parsers_ = [SeqIO.parse(a_file, "fastq") for a_file in files_]
accumulator, denominator = 0, 0

for rec1, rec2, rec3 in zip(parsers_[0], parsers_[1], parsers_[2]):
    check1 = rec1.id == rec3.id
    check2 = rec1.id == rec2.id
    check3 = rec2.id == rec3.id

    if check1 and check2 and check3:
        accumulator+=1
    else:
        print rec1.id, rec2.id, rec3.id
        sys.exit("Headers don't match exactly. Try setting header trim to TRUE")

    denominator+=1

    if denominator == 150:
        print "{} headers were checked and {} were equivalent".format(denominator, accumulator)
        sys.exit(0)

print "{} headers were checked and {} were equivalent".format(denominator, accumulator)



