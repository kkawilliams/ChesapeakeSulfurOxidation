from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sys

to_read = sys.argv[1]
to_write = ""
with open(to_read, 'r') as ifh:
    for idx, line in enumerate(ifh):
        if idx != 0:
            sn, rcbc = line.strip().split("\t")
            bc = str(Seq(rcbc, generic_dna).reverse_complement())
            print "{}: {} -> {}".format(sn, bc, rcbc)
            new_line = "{}\t{}\n".format(sn, bc)
            to_write += new_line
        else:
            pass

with open(to_read+".rc", "w") as ofh:
    ofh.write(to_write)
