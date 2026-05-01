#!/usr/bin/env python3

import sys
from Motifind import MotiFind

if __name__ == "__main__":


	# get user input
	try:
		fasta_filename = sys.argv[1]
		motif_filename = sys.argv[2]
		max_penalty = int(sys.argv[3])
	except IndexError:
		print("Not enough input arguments given.")
		sys.exit(1)
	except ValueError:
		print("please give deviation as integer.")
		sys.exit(1)


	##########################
	##       CASE 1:        ##
	## Scan of targeted seq ##
	##########################

	print("\nCASE 1 - scan of only one sequence in the fasta file")

	# in this example we're only interested in the first fasta entry
	mf1 = MotiFind(fasta_filename, motif_filename, max_penalty)
	target_idx = 0
	mf1.ScanSeqForMotif(target_idx)

	# three ways to proceed with the matches
	print(f"\nresults for entry {target_idx+1} of {fasta_filename}:\n")
	outfile = "results/output_demo.tsv"
	matches = mf1.all_matches	# 1. store them in a variable and use them later in the script
	mf1.ReturnMatches()			# 2. return the results
	mf1.SaveMatches(outfile)	# 3. save them to disk

	input("\n\npress enter to continue to case 2")

	##########################
	##       CASE 2:        ##
	## Fully Automated Scan ##
	##   of the Full FASTA  ##
	##########################

	print("\nCASE 2 - scan a whole fasta file with the possibility to interact with the program")

	# in this example, we want to know all matches in the whole fasta file
	mf2 = MotiFind(fasta_filename, motif_filename, max_penalty)
	mf2.ScanFastaForMotif()

	matches = mf2.all_matches
	mf2.UserInteraction()

	input("\n\npress enter to continue to case 3")

	
	##########################
	##       CASE 3:        ##
	## Extract a Motif from ##
	##   a file for other   ##
	## purposes without the ##
	## need for a .fsa file ##
	##########################

	print("\nCASE 3 - extract a Motif without a FASTA file")

	mf3 = MotiFind(fasta_name=None, motif_name=motif_filename, max_penalty=0)
	motif = mf3.motif

	print("\nExtracted Motif Data Structure:")
	for entry in mf3.motif:
		print(f"  {entry}")

	# if you want to parse the fasta later
	print("\nNow parsing the FASTA file after initialization...")
	mf3.FastaRead(fasta_filename)

	print(f"Successfully loaded {len(mf3.headers)} FASTA headers.")


	print("\n\nEnd of demo. Have a nice day!")





