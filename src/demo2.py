#!/usr/bin/env python3

import sys
from Motifind import MotiFind

if __name__ == "__main__":

	# get user input
	try:
		fasta_filename = sys.argv[1]
		motif_filename = sys.argv[2]
		max_penalty = int(sys.argv[3])
		out_filename_user = sys.argv[4]
	except IndexError:
		print("Not enough input arguments given.")
		sys.exit(1)
	except ValueError:
		print("please give deviation as integer.")
		sys.exit(1)


	mf3 = MotiFind(fasta_filename, motif_filename, max_penalty)
	mf3.ScanFastaForMotif()
	mf3.SaveMatches(out_filename_user)
	# mf3.ReturnMatches()