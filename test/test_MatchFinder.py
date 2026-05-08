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
    mf.sequences = ["GTCACGACTC"]

    # simple motif: A then C
    mf.motif = [
        ('char', {'A'}, 0),
        ('char', {'C'}, 0)
    ]

    # run scan
    mf.ScanSeqForMotif(chosen_idx=0)

    # check result key exists
        assert "test_seq" in mf.all_matches

    # expected matches
    expected = [
        (2, 4, 1.0),
        (6, 8, 1.0)
    ]

    # check number of matches
    assert len(mf.all_matches["test_seq"]) == 2

    # compare exact results
    assert mf.all_matches["test_seq"] == expected


def test_rawseq_simple_match():
    # create object without files
    mf = MotiFind(None, None, max_penalty=5)
    
    # raw sequence string
    rawseq = "ATTTCCAGGA"
    
    # simple motif: C then A
    mf.motif = [
        ('char', {'C'}, 0),
        ('char', {'A'}, 0)
    ]
    
    # run scan
    mf.ScanSeqForMotif(input_sequence=rawseq, head="test_seq")
    
    # check result key exists
    assert "test_seq" in mf.all_matches

    # check exactly one match was found
    assert len(mf.all_matches["test_seq"]) == 1

    # extract match
    start, end, score = mf.all_matches["test_seq"][0]

    # check match position
    assert start == 5
    assert end == 7

    # check confidence score
    assert score == 1.0
        
    