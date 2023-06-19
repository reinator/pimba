#!/usr/bin/env python3
from Bio import SeqIO as seq
from Bio import Phylo as phy
from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
import util
import re

sys.stdout = open( str(snakemake.log), 'w' )
sys.stderr = sys.stdout

def match_and_split( container, search_string ):
    pat = re.compile( search_string.lower() )
    head = set() # aka the matching items
    tail = set() # aka the rest
    for item in container:
        if pat.search( item.lower() ):
            head.add( item )
        else:
            tail.add( item )
    return head, tail

"""
This script takes a constraint tree (newick) containing general search terms at the leafs (such as family names),
and matches them with the taxa in the provided fasta file, such that all matching taxa are attached to the
constraint tree at the matching leaf, generating a multifurcation  of all matches.

Matches are case-insensitive.

Multiple matches are resolved by taking the first that matches.
"""

newick_in   = snakemake.input.tree
fasta_in    = snakemake.input.fasta
newick_out  = snakemake.output[0]

taxa_labels = set()
to_prune    = set()
# first read in the fasta file, and retain all taxa names
for record in seq.parse( fasta_in, "fasta" ):
    taxa_labels.add( record.description )

# next we iterate over all the leaves in the newick tree
tree = list(phy.parse( newick_in, "newick" ))[0]
# for t in tree.get_terminals():
for t in tree.find_clades(terminal=True):
    # and for every leaf we use its label to find matching taxa labels
    matches, taxa_labels = match_and_split( taxa_labels, t.name )
    # if we found any, add them as directly descending leaf nodes
    if len(matches) == 1:
        # special case: if we just have one match, rename the leaf, as we'd otherwise
        # end up with a non-binary tree
        t.name = next(iter( matches ))
    elif len(matches) > 1:
        print(t.name)
        # create new "clades" from the labels and attach to this node
        t.clades = [phy.BaseTree.Clade(name=n) for n in matches]
    else:
        # if we did not find any corresponding taxa, we must prune this leaf,
        # as this would otherwise cause an error with raxml-ng
        to_prune.add( t.name )

for name in to_prune:
    tree.prune( name=name )

# finally, write the result
phy.write(tree, newick_out, "newick", plain=True)
