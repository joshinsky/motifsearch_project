#!/usr/bin/env python3
import os
import sys

######
# Global Functions
######

def usage(msg=None):
    if msg is not None:
        print(msg)
    print("usage: Motifind.py <sequences.fsa> <motif.tsv> <deviation (int)>")
    sys.exit(1)


######
# FastaRead
######

# Reads a fasta file and returns two lists: headers and sequences
def FastaRead(filename):
    # Check if the file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"The file {filename} does not exist.")
    
    headers = []
    sequences = []
    current_seq = []   # Temporary storage for sequence lines

    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()   # Remove whitespace/newlines

                # Skip empty lines
                if not line: 
                    continue

                # Check for header line starting with ">"
                if line.startswith('>'):
                    # Save the in progress sequence, if any
                    if current_seq:
                        sequences.append("".join(current_seq))
                        current_seq = []

                    headers.append(line[1:])   # Store the header without the ">"
                else:
                    # Append sequence line to fragments list
                    current_seq.append(line)

            # Final check to append the last sequence in the file
            if current_seq:
                sequences.append("".join(current_seq))

        # Check for the equality of headers and sequences amount
        if len(headers) != len(sequences):
            raise ValueError("Mismatch between the numbers of headers and sequences.")

        return headers, sequences

    # Error handling for unexpected issues
    except Exception as e:
        raise IOError(f"Error reading FASTA file: {e}")


######
# MotifRead
######

def MotifRead(filename):

    # Check if the file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"The file {filename} does not exist.")

    motif = []
    penalties = []

    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()   # Remove whitespace/newlines

                # Skip empty lines and comments
                if not line or line.startswith('#'): 
                    continue

                parts = line.split('\t')
                if len(parts) < 2:
                    raise IndexError(f"Missing penalty column:\nexpected:\t[letter, penalty]\ngiven:\t{line}")

                motif_letter = parts[0]
                penalty = parts[1]
                
                motif.append(motif_letter)
                penalties.append(penalty)

    # Error handling for unexpected issues
    except IOError as e:
        raise IOError(f"Error reading motif file: {e}")

    # Check for the equality of headers and sequences amount
    if len(motif) != len(penalties):
        raise ValueError("Mismatch between the length of motif and penalties.")

    return motif, penalties

######
# MotifTracker
######





######
# Run Program
######

# get inputs
try:
    fasta_name = sys.argv[1]
    motif_name = sys.argv[2]
    maxpenalty = int(sys.argv[3])
except IndexError:
    usage(msg="Not enough input arguments given.")
except ValueError:
    usage(msg="please give deviation as integer.")

# extract fasta
headers, sequences = FastaRead(fasta_name)
print(headers)
print(sequences)

# extract motif
motif, penalties = MotifRead(motif_name)
# print(motif)
# print(penalties)

