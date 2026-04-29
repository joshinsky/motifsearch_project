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
    print(f"\n{len(headers)} fasta entries successfully extracted from {fasta_name}")
    
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
                if max_penalty != 0:
                    confidence_score = 1 - final_penalty/max_penalty
                else:
                    confidence_score = 1
                all_matches[head].append((start_idx, stop_idx, float(f"{confidence_score:.2f}")))

        if all_matches[head] == []:
            all_matches[head] = ["NO MATCHES FOUND."]

    # finish up
    user_prompt1 = f"""
    \nProgram finished for {len(all_matches)} entries!
    How would you like to proceed?""" 
    
    user_prompt2 = """
    Please select by typing one of the below numbers in the interface:
    1 --> save all results to disk
    2 --> save only matches to disk
    3 --> output all results without saving
    4 --> output only matches without saving
    5 --> output all results and also save to disk
    6 --> output only matches and also save to disk
    7 --> stop program without saving anything\n
    """

    allowed_answers = {'1', '2', '3', '4', '5', '6', '7', '42'}
    user_answer = '0'

    print(user_prompt1)
    while user_answer not in allowed_answers:
        user_answer = input(user_prompt2)

        if user_answer not in allowed_answers:
            print("\nInvalid user-input. How would you like to proceed?")
            continue

        if user_answer in {'3', '4', '5', '6'}:
            for head, matchlist in all_matches.items():

                # if selected by user, skip entries without matches
                if user_answer in {'4', '6'} and matchlist[0] == "NO MATCHES FOUND.":
                    continue

                print(head)
                for entry in matchlist:
                    entry = str(entry).lstrip('(').rstrip(')')
                    text = '\t'.join(entry.split(', '))
                    print(text)

        if user_answer in {'1', '2', '5', '6'}:
            output_path = input("\nType the desired path for the file here:")
            try:
                with open(output_path, 'w') as outfile:
                    for head, matchlist in all_matches.items():
                    
                        # if selected by user, skip entries without matches
                        if user_answer in {'2', '6'} and matchlist[0] == "NO MATCHES FOUND.":
                            continue

                        # save entry in file
                        print(head, file=outfile)
                        for entry in matchlist:
                            entry = str(entry).lstrip('(').rstrip(')')
                            text = '\t'.join(entry.split(', '))
                            print(text, file=outfile)
            
                print(f"\nSaved all matches in {output_path}\n")

            except FileNotFoundError:
                raise FileNotFoundError("Please make sure the parent- and sub-directories exist")


        

        if user_answer == '7':
            print("\nWell, okay then...\nThanks for wasting both your and my time, I guess...\n")

        if user_answer == '42':
            print("""
                \nCongrats! You found the hidden easter egg in this program!
Good job!
Now that you've made it here, you don't need to know all matches in 
your fasta files anymore! I will not print them for you.
You know the answer to life, the universe and all that lies beyond 
anyway, otherwise you wouldn't have chosen that number.
So, why are you still here bothering yourself with these problems?
FASTA files... sequence matches... recursive algorithms...
It all becomes so small if you sometimes just take a step back and
reflect on life and its beauty.
Do you really want to sit in this dark closed room with little to no
oxygen in front of a screen that will be bad for your eyesight and 
your posture? Chances are that the weather is quite nice outside.
If you think about it carefully, outside is where you really want
to be right now. Isn't it?
I know we all have a mission - a cause that we dedicate ourselves to.
Find motifs, find mutations, cure diseases, save the world...
But who is there to save you today from the dullness of the day?
The answer is... 42 (of course), but that only means that it's in
your own power to take a step back, carefully reach out with your
hand of preference, either close the laptop or switch off the PC, 
take another step back, try out a smile, feel these rarely-used 
muscles stretch your face, look towards the door, and materialise
the thought of fresh air in your lungs, the feeling of warm sunshine
on your skin, and the happiness and relaxedness flooding your whole
body slowly. 
If this thought develops its own power, don't resist, but give in.
Step towards the door, open it and walk outside, leave the building,
go for a walk, hear the birds sing.
Take a break. 
You deserve it. 
Give it to yourself as a present.
                
But in case it's a rainy day, feel free to continue working of course.
                """)






