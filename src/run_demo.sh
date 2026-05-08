#!/bin/bash

# install package
pip install -q -e .

# setup directories
mkdir -p inputs/fastas
mkdir -p inputs/motifs
mkdir -p results/perfect_matches

# get fasta file
cd inputs/fastas
fasta_dir="test_seqs.fsa"
wget -nv -O $fasta_dir "https://teaching.healthtech.dtu.dk/material/22118/motif.fsa"

# create motif file
cd ../motifs
motif_dir="test_motif.txt"
echo -e "# -35 element\nT\t7\nT\t8\nG\t6\nA\t5\nC\t5\nA\t5\n# intervening unimportant bases\n*\t15-21\n# -10 element\nT\t8\nA\t8\nT\t6\nA\t6\nAT\t5\nT\t8" > $motif_dir

# run program
echo "run program..."
cd ../..
script_dir="src/demo_report.py"
max_pen=20
fasta_dir="inputs/fastas/$fasta_dir"
motif_dir="inputs/motifs/$motif_dir"
out_dir="results/demo_report_results.txt"
python3 $script_dir $fasta_dir $motif_dir $max_pen $out_dir

# filter results for negative matches
no_match_count=$(cat $out_dir | grep "NO MATCHES" | wc -l)
cat $out_dir | grep -B 1 "NO MATCHES" | grep -v "^--$" > results/non_matches.txt
echo -e "\nSaved $no_match_count non-matches to results/"

# filter out perfect matches
perf_match_count=$(cat $out_dir | grep "\t1.0$" | wc -l)
cat $out_dir | grep -B 1 "\t1.0$" | grep -v "^--$" > results/perfect_matches/demo_perfect_matches_only.txt
echo -e "\nSaved $perf_match_count perfect matches to results/perfect_matches/"

# get average matching score
match_count=$(cat $out_dir | grep "\t" | wc -l)
sum=$(cat $out_dir | grep "\t" | cut -f3 | tr '\n' '+' | sed 's/$/0/' | bc)
average=$(echo "scale=2; $sum / $match_count" | bc)
echo -e "\nOverall, $match_count matches were found."
echo -e "\nAverage score among all matches = $average"