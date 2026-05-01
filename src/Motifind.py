#!/usr/bin/env python3
import os
import sys


class MotiFind:

    def __init__(self, fasta_name, motif_name, max_penalty):
        self.fasta_name = fasta_name
        self.motif_name = motif_name
        self.max_penalty = max_penalty

        self.headers = []
        self.sequences = []
        self.motif = []
        self.all_matches = {}


    ######
    # Global Functions
    ######

    def usage(self, msg=None):
        if msg is not None:
            print(msg)
        print("usage: Motifind.py <sequences.fsa> <motif.tsv> <deviation (int)>")
        sys.exit(1)


    ###############
    ## FastaRead ##
    ###############

    # Reads a fasta file and returns two lists: headers and sequences
    def FastaRead(self, filename):

        # Check if the file exists
        if not os.path.exists(filename):
            raise FileNotFoundError(f"The file {filename} does not exist.")
        
        headers = []
        sequences = []
        current_seq = []   # Temporary storage for sequence lines

        try:
            with open(filename, 'r') as f:
                for line in f:
                    line = line.strip()     # Remove whitespace/newlines

                    # Skip empty lines
                    if not line: 
                        continue

                    # Check for header line starting with ">"
                    if line.startswith('>'):
                        # Save the in progress sequence, if any, then save the header
                        if current_seq:
                            sequences.append("".join(current_seq))
                            current_seq = []
                        headers.append(line[1:])   # Exclude '>' from headers

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

    def MotifRead(self, filename):

        # Check if the file exists
        if not os.path.exists(filename):
            raise FileNotFoundError(f"The file {filename} does not exist.")

        motif = []

        try:
            with open(filename, 'r') as f:
                for i, line in enumerate(f):
                    line = line.strip()     # Remove whitespace/newlines

                    # Skip empty lines and comments
                    if not line or line.startswith('#'): 
                        continue

                    # extract columns of the .tsv
                    parts = line.split('\t')
                    if len(parts) < 2:
                        raise IndexError(f"Error in line {i} of {filename}: Missing second column:\nexpected:\t[letter, penalty]\ngiven:\t{line}")

                    # Gap case
                    if parts[0] == '*':
                        entry_type = 'gap'

                        # Allow both '-' and '_' as separators
                        parts[1] = parts[1].replace('_', '-')

                        if '-' in parts[1]:
                            gap_bounds = parts[1].split('-')

                            if len(gap_bounds) != 2:
                                raise ValueError(f"Error in line {i} of {filename}: gap must be in format 'min-max'. Got: {parts[1]}")

                            # convert bounds to int
                            try:
                                lower_bound = int(gap_bounds[0])
                                upper_bound = int(gap_bounds[1])
                            except ValueError:
                                raise ValueError(f"Error in line {i} of {filename}: gap bounds should be given as integers.\n{gap_bounds} was given which are not integers.")

                            if lower_bound > upper_bound:
                                raise ValueError(f"Error in line {i} of {filename}: lower bound should not be > upper bound in gap: {parts[1]}")

                        # case: gap spans only one single value
                        else:
                            try:
                                lower_bound = upper_bound = int(parts[1])
                            except ValueError:
                                raise ValueError(f"Error in line {i} of {filename}: gap value must be an integer. Got: {parts[1]}")

                        motif_entry = (entry_type, lower_bound, upper_bound)
                        
                    # Character case
                    else:
                        entry_type = 'char'
                        motif_chars = set(parts[0])
                        penalty = parts[1]

                        # convert penalty to float
                        try:
                            penalty = float(penalty)
                        except ValueError:
                            raise ValueError(f"Error in line {i} of {filename}: penalties should be given as numbers.\nYour input '{penalty}' is non-numerical.")

                        if penalty < 0:
                            raise ValueError(f"Error in line {i} of {filename}: penalties should not be negative numbers. {penalty} is negative.")

                        motif_entry = (entry_type, motif_chars, penalty)
                    
                    motif.append(motif_entry)

        # Error handling for unexpected issues
        except IOError as e:
            raise IOError(f"Error reading motif file: {e}")

        return motif


    #################
    ## MatchFinder ##
    #################

    def MatchFinder(self, seq, current_seq_idx, motif, current_motif_idx, accum_penalty, max_penalty):

        # case: we exceded the end of the seq without a match
        if current_seq_idx > len(seq):
            matchlist = []
            return matchlist

        # case: we reached the end of the motif without crossing the max_penalty threshold --> match!
        if current_motif_idx >= len(motif):
            matchlist = [[current_seq_idx, accum_penalty]]
            return matchlist

        # case: we reached the end of the motif without a match
        if current_seq_idx == len(seq):
            matchlist = []
            return matchlist

        # case: our next motif element is a defined character
        if motif[current_motif_idx][0] == 'char':
            seq_char = seq[current_seq_idx]
            motif_chars = motif[current_motif_idx][1]
            mismatch_penalty = motif[current_motif_idx][2]

            # in a mismatch, increase the total penalty
            if seq_char not in motif_chars:
                accum_penalty += mismatch_penalty

                # if our maximum penalty was exceded, we return "no match"
                if accum_penalty > max_penalty:
                    matchlist = []
                    return matchlist

            # recursive function call on the next motif and sequence position
            current_seq_idx += 1
            current_motif_idx += 1  
            return self.MatchFinder(seq, current_seq_idx, motif, current_motif_idx, accum_penalty, max_penalty)

        # case: our next motif elements can be any character
        elif motif[current_motif_idx][0] == 'gap':
            min_gap_len = motif[current_motif_idx][1]
            max_gap_len = motif[current_motif_idx][2] + 1
            
            current_motif_idx += 1
            gap_matches = []

            # recursively try out all different gap scenarios and collect any occuring matches
            for i in range(min_gap_len, max_gap_len):
                gap_accum_penalty = accum_penalty
                gap_seq_idx = current_seq_idx + i
                matchlist = self.MatchFinder(seq, gap_seq_idx, motif, current_motif_idx, gap_accum_penalty, max_penalty)
                gap_matches.extend(matchlist[:])

            return gap_matches


    ##################
    ## Main Program ##
    ##################

    def run(self):

        # get user inputs
        fasta_name = self.fasta_name
        motif_name = self.motif_name
        max_penalty = self.max_penalty

        # extract fasta
        headers, sequences = self.FastaRead(fasta_name)
        print(f"\n{len(headers)} fasta entries successfully extracted from {fasta_name}")
        
        # extract motif
        motif = self.MotifRead(motif_name)
        print(f"The motif provided in {motif_name} has a minimum length of {len(motif)}.")
        print(f"So, let's get started with a maximum mismatch score of {max_penalty}!")

        # check out all given fasta entries, one after another
        all_matches = {}
        for i in range(len(sequences)):
            seq = sequences[i]
            head = headers[i]
            all_matches[head] = []
            
            # loop over all possible sequence start indeces and search for motif matches
            start_idx = 0
            for start_idx in range(len(seq)):
                motif_idx = 0
                start_penalty = 0
                matches = self.MatchFinder(seq, start_idx, motif, motif_idx, start_penalty, max_penalty)
            
                # store any matches found in a dict {head: [(startidx1, stopidx1, confidence1), ...]}
                for match in matches:
                    stop_idx = match[0]
                    final_penalty = match[1]
                    if max_penalty != 0:
                        confidence_score = 1 - final_penalty/max_penalty
                    else:
                        confidence_score = 1
                    all_matches[head].append((start_idx, stop_idx, float(f"{confidence_score:.2f}")))

            # take note of each entry that didn't yield matches
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
        user_answer = None

        print(user_prompt1)
        while user_answer not in allowed_answers:
            user_answer = input(user_prompt2)

            if user_answer not in allowed_answers:
                print("\nInvalid user-input. How would you like to proceed?")
                continue

            # case: user wants to print the output on screen
            if user_answer in {'3', '4', '5', '6'}:
                for head, matchlist in all_matches.items():

                    if user_answer in {'4', '6'} and matchlist[0] == "NO MATCHES FOUND.":
                        continue

                    print(head)
                    for entry in matchlist:
                        entry = str(entry).lstrip('(').rstrip(')')
                        text = '\t'.join(entry.split(', '))
                        print(text)

            # case: user wants to store the ouput on disk
            if user_answer in {'1', '2', '5', '6'}:
                output_path = input("\nType the desired path for the file here:")
                try:
                    with open(output_path, 'w') as outfile:
                        for head, matchlist in all_matches.items():

                            if user_answer in {'2', '6'} and matchlist[0] == "NO MATCHES FOUND.":
                                continue

                            print(head, file=outfile)
                            for entry in matchlist:
                                entry = str(entry).lstrip('(').rstrip(')')
                                text = '\t'.join(entry.split(', '))
                                print(text, file=outfile)

                    print(f"\nSaved all matches in {output_path}\n")

                except FileNotFoundError:
                    raise FileNotFoundError("Please make sure the parent- and sub-directories exist")

            # case: user doesn't want to do anything
            if user_answer == '7':
                print("\nWell, okay then...\nThanks for wasting both your and my time, I guess...\n")

            # case: user finds the easter egg
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

                rainy_day = input(f"\nIs it a rainy day? (y/n)")
                if rainy_day.lower() == 'y':
                    user_answer = None
                    continue
                elif rainy_day.lower() == 'n':
                    print(f"That's great to hear! Enjoy the day!")
                    sys.exit(1)
                else:
                    print(f"I don't understand that answer. That must mean you're very overworked and should take a break!")
                    sys.exit(1)


# ---- MAIN ----
if __name__ == "__main__":

    try:
        fasta_name = sys.argv[1]
        motif_name = sys.argv[2]
        max_penalty = int(sys.argv[3])
    except IndexError:
        print("Not enough input arguments given.")
        sys.exit(1)
    except ValueError:
        print("please give deviation as integer.")
        sys.exit(1)

    program = MotiFind(fasta_name, motif_name, max_penalty)
    program.run()