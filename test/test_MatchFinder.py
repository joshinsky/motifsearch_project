import pytest
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), "../src"))

from Motifind import MotiFind

def test_scanseq_simple_match():
    # create object without files 
    mf = MotiFind(None, None, max_penalty=5)

    # manually defined FASTA-like data
    mf.headers = ["test_seq"]
    mf.sequences = ["AC"]

    # simple motif: A then C
    mf.motif = [
        ('char', {'A'}, 0),
        ('char', {'C'}, 0)
    ]

    # run scan
    mf.ScanSeqForMotif(chosen_idx=0)

    # check results exist
    assert "test_seq" in mf.all_matches

    # check at least one match found
    assert len(mf.all_matches["test_seq"]) > 0

    # check confidence is perfect (no penalty)
    start, end, score = mf.all_matches["test_seq"][0]
    assert score == 1.0
    assert start == 0
    assert end == 2
