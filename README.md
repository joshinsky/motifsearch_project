# MotiFind

**MotiFind** is a robust, rule-based exact motif search package for bioinformatics. 

Discovering conserved sequence patterns (motifs) is a computationally difficult task, especially when dealing with the low signal-to-noise ratio in biological datasets where motifs are often "planted" within large regulatory regions. While many tools rely on heuristic or probabilistic models (like HMMs or PSSMs), MotiFind utilizes a deterministic, recursive backtracking algorithm. 

It allows users to explicitly define allowed characters, assign custom mismatch penalties, and incorporate highly flexible gap regions to model variable spacing. This makes it an ideal, transparent tool for identifying transcription factor binding sites (TFBS), diagnostic probes, and protein families across massive FASTA datasets.

## Key Features
* **Custom Mismatch Penalties:** Define exact penalty scores for specific position mismatches.
* **Variable Gap Handling:** Recursively explores dynamic gap length ranges (e.g., 15-21 intervening bases).
* **Object-Oriented & Modular:** Easily importable into larger bioinformatics pipelines or Jupyter Notebooks.
* **Interactive CLI:** Built-in methods for cleanly outputting, filtering, and saving matches to `.txt` or `.tsv` files.

---

## Quick Start Guide

### 1. Installation
MotiFind is structured as a standard Python package. It is recommended to install it in editable mode within your virtual environment so all internal module dependencies resolve correctly.

```bash
# Clone the repository
git clone https://github.com/joshinsky/motifsearch_project.git
cd motifsearch_project

# Install the package in editable mode
pip install -e .
```

### 2. Running the Program
You can run the built-in demonstration script from your terminal to see MotiFind in action immediately. The demo automatically downloads a sample FASTA file, generates a motif definition file, and runs the search:

```bash
# Run the automated bash pipeline
./src/demo/run_demo.sh
```

Alternatively, you can run the core Python script directly by passing your own FASTA file, motif definition file, maximum penalty threshold, and desired output path:

```bash
# Direct python execution
mkdir results
python3 src/demo/demo2.py inputs/fastas/test_seqs.fsa inputs/motifs/test_motif.txt 15 results/my_results.txt
```

### 3. Usage in Python
Because it is installed as a package, you can easily import and use ** **MotiFind** in any custom Python script:

```Python
from motifind import MotiFind

# Initialize with files and a max penalty threshold of 5
mf = MotiFind(
    fasta_name="data.fsa",
    motif_name="signal.txt",
    max_penalty=5
)

# Scan the entire FASTA dataset
mf.ScanFastaForMotif()

# Save the results
mf.SaveMatches("results.tsv")
``
