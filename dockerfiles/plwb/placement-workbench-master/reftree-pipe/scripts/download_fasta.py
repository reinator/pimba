import csv, sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
import util
import pandas as pd

from Bio import Entrez as gb
from Bio import SeqIO as seq
gb.email = "pierre@barbera-bio.info"

import re
genbank_acc_pattern 		= re.compile("^([A-Z]{1,6}_?[0-9]{2,8}.?[0-9]*)$")
genbank_acc_pattern_version = re.compile("^([A-Z]{1,6}_?[0-9]{2,8}\.[0-9]+)$")

def dl_to_fasta( accession_file, fasta_dest_path ):

	# read in file containing acession numbers / labels
	accessions = pd.read_table(accession_file, dtype=str, sep=',')

	label_col = snakemake.config["data"]["accession_file_structure"]["label_column"]
	acc_col = snakemake.config["data"]["accession_file_structure"]["accession_column"]

	# label_col = "label"
	# acc_col = "accession"

	# validate the accessions
	for i, row in accessions.iterrows():
		acc = row[ acc_col ]
		if( not genbank_acc_pattern.match(acc) ):
				util.fail( "'{acc}' does not look like a valid accession" )
		elif( not genbank_acc_pattern_version.match(acc) ):
			util.warn( "'{acc}' does not seem to have a version, appending '.1'" )
			accessions.at[ i, acc_col ] = acc + ".1"

	# now we index such that we can better associate accessions with their labels
	accessions.set_index( acc_col, inplace=True )

	# request accessions from genbank
	try:
		gb_handle = gb.efetch(
			db="nucleotide",
			id=accessions.index.tolist(),
			rettype="fasta",
			retmode="text")
		with open( fasta_dest_path, 'w+' ) as out_fasta:
			# write out to (compressed?) fasta
			for record in seq.parse( gb_handle, "fasta" ):
				# replace the ID with the label specified in the input csv
				label = accessions.at[ record.id, label_col ]
				record.id = label
				record.description = ""
				out_fasta.write( record.format("fasta") )

		gb_handle.close()
	except IOError:
		util.fail("Could not connect to GenBank.")

def copy_tsv( input_csv_path, output_fasta_path ):
	input_tsv_path = "{}/{}.tsv".format( util.dirname(input_csv_path), util.filename(input_csv_path) )
	output_tsv_path = "{}/{}.tsv".format( util.dirname(output_fasta_path), util.filename(output_fasta_path) )

	util.expect_file_exists( input_tsv_path )

	util.copy( input_tsv_path, output_tsv_path )

if __name__ == "__main__":
	with open(snakemake.log[0], "w") as f:
		sys.stderr = sys.stdout = f
		# first try to copy over the taxonomy tsv files, so we can short-circuit if that fails
		# if snakemake.params.copy_tsvs_over:
		# 	copy_tsv( snakemake.input[0], snakemake.output[0] )

		dl_to_fasta( snakemake.input[0], snakemake.output[0] )
