#!/usr/bin/env python3
import os
import sys

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






######
# MotifTracker
######


