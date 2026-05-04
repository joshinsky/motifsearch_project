#!/usr/bin/env python3
import os
import sys


class MotiFind:

    def __init__(self, fasta_name, motif_name, max_penalty):
        self.fasta_name = fasta_name
        self.motif_name = motif_name
        self.max_penalty = max_penalty

        self.sequences = []
        self.headers = []
        self.motif = []
        self.all_matches = {}

        # extract headers, sequences and motif if available
        if fasta_name is not None:
            self.FastaRead()
        if motif_name is not None:
            self.MotifRead()

        self.all_matches = {}


    ###############
    ## FastaRead ##
    ###############

    # Reads a fasta file and returns two lists: headers and sequences
    def FastaRead(self, fasta_filename=None):

        if fasta_filename is None:
            filename = self.fasta_name
        else:
            filename = fasta_filename


        if filename is None:
            raise ValueError("No FASTA filename provided.")

        # Check if the file exists
        if not os.path.exists(filename):
            raise FileNotFoundError(f"The file {filename} does not exist.")
        
        self.headers = []
        self.sequences = [] 

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
                            self.sequences.append("".join(current_seq))
                            current_seq = []
                        self.headers.append(line[1:])   # Exclude '>' from headers

                    else:
                        # Append sequence line to fragments list
                        current_seq.append(line)

                # Final check to append the last sequence in the file
                if current_seq:
                    self.sequences.append("".join(current_seq))

            # Check for the equality of headers and sequences amount
            if len(self.headers) != len(self.sequences):
                raise ValueError("Mismatch between the numbers of headers and sequences.")

            print(f"\n{len(self.headers)} fasta entries successfully extracted from {filename}\n")
            return self.headers, self.sequences

        # Error handling for unexpected issues
        except OSError as e:
            raise IOError(f"Error reading FASTA file: {e}")


    ###############
    ## MotifRead ##
    ###############

    def MotifRead(self, motif_filename=None):

        if motif_filename == None:
            filename = self.motif_name
        else:
            filename = motif_filename

        # Check if the file exists
        if not os.path.exists(filename):
            raise FileNotFoundError(f"The file {filename} does not exist.")

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
                        raise IndexError(f"\n\nERROR IN LINE {i} of {filename}: Missing second column:\nexpected:\t[letter, penalty]\ngiven:\t{line}\n")

                    # Gap case
                    if parts[0] == '*':
                        entry_type = 'gap'

                        # Allow both '-' and '_' as separators
                        parts[1] = parts[1].replace('_', '-')

                        if '-' in parts[1]:
                            gap_bounds = parts[1].split('-')

                            if len(gap_bounds) != 2:
                                raise ValueError(f"\n\nERROR IN LINE {i} of {filename}: gap must be in format 'min-max'. Got: {parts[1]}\n")

                            # convert bounds to int
                            try:
                                lower_bound = int(gap_bounds[0])
                                upper_bound = int(gap_bounds[1])
                            except ValueError:
                                raise ValueError(f"\n\nERROR IN LINE {i} of {filename}: gap bounds should be given as integers.\n{gap_bounds} was given which are not integers.\n")

                            if lower_bound > upper_bound:
                                raise ValueError(f"\n\nERROR IN LINE {i} of {filename}: lower bound should not be > upper bound in gap: {parts[1]}\n")

                        # case: gap spans only one single value
                        else:
                            try:
                                lower_bound = upper_bound = int(parts[1])
                            except ValueError:
                                raise ValueError(f"\n\nERROR IN LINE {i} of {filename}: gap value must be an integer. Got: {parts[1]}\n")

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
                            raise ValueError(f"\n\nERROR IN LINE {i} of {filename}: penalties should be given as numbers.\nYour input '{penalty}' is non-numerical.\n")

                        if penalty < 0:
                            raise ValueError(f"\n\nERROR IN LINE {i} of {filename}: penalties should not be negative numbers. {penalty} is negative.\n")

                        motif_entry = (entry_type, motif_chars, penalty)
                    
                    self.motif.append(motif_entry)

            return self.motif

        # Error handling for unexpected issues
        except IOError as e:
            raise IOError(f"Error reading motif file: {e}")


    #################
    ## MatchFinder ##
    #################

    def MatchFinder(self, seq, current_seq_idx, current_motif_idx, accum_penalty):

        max_penalty = self.max_penalty
        motif = self.motif

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
            return self.MatchFinder(seq, current_seq_idx, current_motif_idx, accum_penalty)

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
                matchlist = self.MatchFinder(seq, gap_seq_idx, current_motif_idx, gap_accum_penalty)
                gap_matches.extend(matchlist[:])

            return gap_matches


    ###################
    ## Sequence Scan ##
    ###################

    def ScanSeqForMotif(self, chosen_idx:int = None, input_sequence:str = None, head:str = None):

        # Enable inputting a sequence or an index
        if input_sequence is not None and chosen_idx is not None:
            raise ValueError("Provide either input_sequence or chosen_idx, not both.")

        if input_sequence is None and chosen_idx is None:
            raise ValueError("You must provide either input_sequence or chosen_idx")

        if input_sequence is not None:
            if not isinstance(input_sequence, str):
                raise TypeError("input_sequence must be a string")

            seq = input_sequence
            head = head or "input_sequence"

        else:
            if not isinstance(chosen_idx, int):
                raise TypeError("chosen_idx must be an integer")

            try:
                seq = self.sequences[chosen_idx]
                head = self.headers[chosen_idx]

            except IndexError:
                raise IndexError("chosen_idx is out of range.")

        self.all_matches[head] = []

        # loop over all possible sequence start indeces and search for motif matches
        start_idx = 0
        for start_idx in range(len(seq)):
            motif_idx = 0
            start_penalty = 0
            matches = self.MatchFinder(seq, start_idx, motif_idx, start_penalty)

            # store any matches found in a dict {head: [(startidx1, stopidx1, confidence1), ...]}
            for match in matches:
                stop_idx = match[0]
                final_penalty = match[1]
                if self.max_penalty != 0:
                    confidence_score = 1 - final_penalty/self.max_penalty
                else:
                    confidence_score = 1
                self.all_matches[head].append((start_idx, stop_idx, float(f"{confidence_score:.2f}")))

        # take note of each entry that didn't yield matches
        if not self.all_matches[head]:
            self.all_matches[head] = ["NO MATCHES FOUND."]

    def ScanFastaForMotif(self):

        # check out all given fasta entries, one after another
        for i in range(len(self.sequences)):
            self.ScanSeqForMotif(i)

        return self.all_matches


    ############################
    ## Save or Return Matches ##
    ############################

    def SaveMatches(self, output_path, matches_only=False):
        try:
            with open(output_path, 'w') as outfile:
                for head, matchlist in self.all_matches.items():
                    if matches_only and matchlist[0] == "NO MATCHES FOUND.":
                        continue

                    print(head, file=outfile)
                    for entry in matchlist:
                        entry = str(entry).lstrip('(').rstrip(')')
                        text = '\t'.join(entry.split(', '))
                        print(text, file=outfile)

            print(f"\nSaved all matches in {output_path}\n")

        except FileNotFoundError:
            raise FileNotFoundError("Please make sure the parent- and sub-directories exist")

    def ReturnMatches(self, matches_only=False):
        for head, matchlist in self.all_matches.items():
            if matches_only and matchlist[0] == "NO MATCHES FOUND.":
                continue
            print(head)
            for entry in matchlist:
                entry = str(entry).lstrip('(').rstrip(')')
                text = '\t'.join(entry.split(', '))
                print(text)

    def UserInteraction(self):
        user_prompt1 = f"""
        \nProgram finished! How would you like to proceed?""" 
        
        user_prompt2 = """
        Please select by typing one of the below numbers in the interface:
        1 --> save all results to disk
        2 --> save only matches to disk
        3 --> output all results without saving
        4 --> output only matches without saving
        5 --> output all results and also save to disk
        6 --> output only matches and also save to disk
        7 --> stop program\n
        """

        user_prompt3 = "\nWould you like to do anything else?"

        print(user_prompt1)

        user_answer = None
        first_time = True
        while user_answer != '7':
            user_answer = input(user_prompt2)

            if user_answer not in {'1', '2', '3', '4', '5', '6', '7', '42'}:
                print("\nInvalid user-input. How would you like to proceed?")
                continue

            # does the user only want matches or all results
            if user_answer in {'2', '4', '6'}:
                only_matches = True
            else:
                only_matches = False

            # case: user wants to print the output on screen
            if user_answer in {'3', '4', '5', '6'}:
                self.ReturnMatches(matches_only=only_matches)
                first_time = False

            # case: user wants to store the ouput on disk
            if user_answer in {'1', '2', '5', '6'}:
                output_path = input("\nType the desired path for the file here:")
                self.SaveMatches(output_path, matches_only=only_matches)
                first_time = False

            # case: user doesn't want to do anything
            if user_answer == '7':
                if first_time:
                    print("\nWell, okay then...\nThanks for wasting both your and my time, I guess...\n")
                else:
                    print("\nAlrighty!\nThanks for stopping by!\n")
                continue

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
                    print("\nOkay then! Back to work!")
                    user_answer = None
                    continue
                elif rainy_day.lower() == 'n':
                    print(f"\nThat's great to hear! Enjoy the day!")
                    sys.exit(1)
                else:
                    print(f"\nI don't understand that answer. That must mean you're very overworked and should take a break!")
                    sys.exit(1)

            print(user_prompt3)
