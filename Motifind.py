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


###############
## FastaRead ##
###############

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
                    # Save the in progress sequence, if any, then save the header
                    if current_seq:
                        sequences.append("".join(current_seq))
                        current_seq = []
                    headers.append(line)

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



###############
## MotifRead ##
###############

def MotifRead(filename):

    # Check if the file exists
    if not os.path.exists(filename):
        raise FileNotFoundError(f"The file {filename} does not exist.")

    motif = []

    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()   # Remove whitespace/newlines

                # Skip empty lines and comments
                if not line or line.startswith('#'): 
                    continue

                parts = line.split('\t')
                if len(parts) < 2:
                    raise IndexError(f"Missing second column:\nexpected:\t[letter, penalty]\ngiven:\t{line}")

                if parts[0] == '*':
                    entry_type = 'gap'
                    gap_bounds = parts[1].split('-')
                    try:
                        lower_bound = int(gap_bounds[0])
                        upper_bound = int(gap_bounds[1])
                    except ValueError:
                        raise ValueError(f"Error: gap bounds should be given as integers.\n{gap_bounds} was given which are not integers.")    
                    motif_entry = (entry_type, lower_bound, upper_bound)
                    
                else:
                    entry_type = 'char'
                    motif_chars = set(parts[0])
                    penalty = parts[1]
                    try:
                        penalty = float(penalty)
                    except ValueError:
                        raise ValueError(f"Error: penalties should be given as numbers.\nYour input '{penalty}' is non-numerical.")
                    motif_entry = (entry_type, motif_chars, penalty)
                
                motif.append(motif_entry)

    # Error handling for unexpected issues
    except IOError as e:
        raise IOError(f"Error reading motif file: {e}")

    return motif



#################
## MatchFinder ##
#################

def MatchFinder(seq, current_seq_idx, motif, current_motif_idx, accum_penalty, max_penalty):

    if current_seq_idx > len(seq):
        matchlist = []
        return matchlist

    if current_motif_idx >= len(motif):
        matchlist = [[current_seq_idx, accum_penalty]]
        return matchlist

    if current_seq_idx == len(seq):
        matchlist = []
        return matchlist

    if motif[current_motif_idx][0] == 'char':
        seq_char = seq[current_seq_idx]
        motif_chars = motif[current_motif_idx][1]
        mismatch_penalty = motif[current_motif_idx][2]

        if seq_char not in motif_chars:
            accum_penalty += mismatch_penalty
            if accum_penalty > max_penalty:
                matchlist = []
                return matchlist

        current_seq_idx += 1
        current_motif_idx += 1  
        return MatchFinder(seq, current_seq_idx, motif, current_motif_idx, accum_penalty, max_penalty)

    elif motif[current_motif_idx][0] == 'gap':
        min_gap_len = motif[current_motif_idx][1]
        max_gap_len = motif[current_motif_idx][2] + 1
        
        current_motif_idx += 1
        gap_matches = []
        for i in range(min_gap_len, max_gap_len):
            gap_accum_penalty = accum_penalty
            gap_seq_idx = current_seq_idx+i
            matchlist = MatchFinder(seq, gap_seq_idx, motif, current_motif_idx, gap_accum_penalty, max_penalty)
            gap_matches.extend(matchlist[:])
        return gap_matches


##################
## Main Program ##
##################

if __name__ == "__main__":

    # get user inputs
    try:
        fasta_name = sys.argv[1]
        motif_name = sys.argv[2]
        max_penalty = int(sys.argv[3])
    except IndexError:
        usage(msg="Not enough input arguments given.")
    except ValueError:
        usage(msg="please give deviation as integer.")

    # extract fasta
    headers, sequences = FastaRead(fasta_name)
    print(f"{len(headers)} fasta entries successfully extracted from {fasta_name}")
    
    # extract motif
    motif = MotifRead(motif_name)
    print(f"The motif provided in {motif_name} has a minimum length of {len(motif)}.")
    print(f"So, let's get started with a maximum mismatch score of {max_penalty}!")

    # find matches
    all_matches = {}
    for i in range(len(sequences)):
        seq = sequences[i]
        head = headers[i]
        all_matches[head] = []
        
        start_idx = 0
        for start_idx in range(len(seq)):
            motif_idx = 0
            start_penalty = 0
            matches = MatchFinder(seq, start_idx, motif, motif_idx, start_penalty, max_penalty)
        
            for match in matches:
                stop_idx = match[0]
                final_penalty = match[1]
                all_matches[head].append((start_idx, stop_idx, final_penalty))

    # finish up
    user_prompt1 = f"""
    \nProgram finished!
    \nFound one or more matches for {len(all_matches)} entries!
    \nHow would you like to proceed?
    """ 
    
    user_prompt2 = """
    \nPlease select by typing one of the below numbers in the interface:
    \n1 --> save matches to disk
    \n2 --> output matches without saving
    \n3 --> output matches and save to disk
    \n4 --> stop program without saving matches
    """
    allowed_answers = {'1', '2', '3', '4', '42'}
    user_answer = '0'

    print(user_prompt1)
    while user_answer not in allowed_answers:
        user_answer = input(user_prompt2)

        if user_answer in {'1', '3'}:
            output_filename = input("Type the desired name for the file here (without file-suffix):")
            with open(f"{output_filename}.tsv", 'w') as outfile:
                for head, matchlist in all_matches.items():
                    print(head, file=outfile)
                    for match in matchlist:
                        print(match, file=outfile)
        
        if user_answer in {'2', '3'}:
            for head, matchlist in all_matches.items():
                print(head)
                for match in matchlist:
                    print(match)

        if user_answer == '4':
            print("\nWell, okay then...\nThanks for wasting both your and my time, I guess...")

        if user_answer == '42':
            print("""
                \nCongrats! You found the hidden easter egg in this program!
                \nGood job!
                \nNow that you've made it here, you don't need to know all matches in 
                \nyour fasta files anymore! I will not print them for you.
                \nYou know the answer to life, the universe and all that lies beyond 
                \nanyway, otherwise you wouldn't have chosen this number.
                \nSo, why are you still here bothering yourself with these problems?
                \nFASTA files... sequence matches... recursive algorithms...
                \nIt all becomes so small if you sometimes just take a step back and
                \nreflect on life and its beauty.
                \nDo you really want to sit in this dark closed room with little to no
                \noxygen in front of a screen that will be bad for your eyesight and 
                \nyour posture? Chances are that the weather is quite nice outside.
                \nIf you think about it carefully, outside is where you really want
                \nto be right now. Isn't it?
                \nI know we all have a mission - a cause that we dedicate ourselves to.
                \nFind motifs, find mutations, cure diseases, save the world...
                \nBut who is there to save you today from the dullness of the day?
                \nThe answer is... 42 (of course), but that only means that it's in
                \nyour own power to take a step back, carefully reach out with your
                \nhand of preference, either close the laptop or switch off the PC, 
                \ntake another step back, try out a smile, feel these rarely-used 
                \nmuscles stretch your face, look towards the door, and materialise
                \nthe thought of fresh air in your lungs, the feeling of warm sunshine
                \non your skin, and the happiness and relaxedness flooding your whole
                \nbody slowly. 
                \nIf this thought develops its own power, don't resist, but give in.
                \nStep towards the door, open it and step outside, leave the building,
                \ngo for a walk, hear the birds sing.
                \nTake a break. You deserve it. Give it to yourself as a present.
                \n
                \nBut in case it's a rainy day, you can feel free to continue working.
                """)
            sys.exit(0)






