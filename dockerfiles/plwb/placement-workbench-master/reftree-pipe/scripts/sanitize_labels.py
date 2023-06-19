#!/usr/bin/env python3

from Bio import SeqIO as seq
from Bio import Phylo as phy
import sys
import re

in_file = sys.argv[1]
out_file = sys.argv[2]
mode = sys.argv[3]

# sanitize the extra complicated via utf8, converting each weird opcode to its own underscore
# because thats how IQTREE does it, yay
def sanitize( label ):
    rec = label.encode('utf-8')
    rec = re.sub(b'[^\w]|\\[\w]{3}', b'_', rec)
    rec = rec.decode('utf-8')
    if rec != label:
        print(f"{label} => {rec}")
    return rec

if __name__ == "__main__":
    if mode == "fasta":
        with open(out_file, 'w') as of:
            for record in seq.parse(in_file, "fasta"):
                record.id = sanitize( record.id )
                record.description = record.id
                seq.write(record, of, "fasta")
    elif mode == "newick":
        trees = list(phy.parse(in_file, "newick"))
        for tree in trees:
            for t in tree.get_terminals():
                t.name = sanitize(t.name)
        phy.write(trees, out_file, "newick")
    else:
        raise Exception( "Unsupported mode: " + mode )
