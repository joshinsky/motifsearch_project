#!/bin/bash

# install package
pip install -q -e .

# setup directories
mkdir -p inputs/fastas
mkdir -p inputs/motifs
mkdir -p results/perfect_matches

# get fasta file
cd inputs/fastas
wget -nv -O test_seqs.fsa "https://teaching.healthtech.dtu.dk/material/22118/motif.fsa"

# create motif file
cd ../motifs
echo -e "# -35 element\nT\t7\nT\t8\nG\t6\nA\t5\nC\t5\nA\t5\n# intervening unimportant bases\n*\t15-21\n# -10 element\nT\t8\nA\t8\nT\t6\nA\t6\nAT\t5\nT\t8" #> test_motif.txt

# run program
cd ../..
max_pen=15
fasta_dir="inputs/fastas/test_seqs.fsa"
motif_dir="inputs/motifs/test_motif.txt"
out_dir="results/demo_results.txt"
python3 src/demo2.py $fasta_dir $motif_dir $max_pen $out_dir

# filter results for negative matches
echo -e "\n\nLet's filter out only non-matches from $out_dir"
cat $out_dir | grep -B 1 "NO MATCHES" | grep -v "^--$" > results/non_matches.txt

# filter out perfect matches
echo -e "\nLet's filter out only perfect matches"
cat $out_dir | grep -B 1 "\t1.0$" | grep -v "^--$" > results/perfect_matches/demo_perfect_matches_only.txt

# get average matching score
sum=$(cat $out_dir | grep "\t" | cut -f3 | tr '\n' '+' | sed 's/$/0/' | bc)
count=$(cat $out_dir | grep "\t" | wc -l)
average=$(echo "scale=2; $sum / $count" | bc)
echo -e "Average score among all matches = $average"