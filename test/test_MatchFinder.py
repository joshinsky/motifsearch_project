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
    assert set(mf.all_matches["test_seq"]) == set(expected)


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

    # check match position and confidence score
    assert (start, end, score) == (5, 7, 1.0)
    

def test_fasta_missing_file():
    with pytest.raises(FileNotFoundError):
        MotiFind("doesnt_exist.fsa", None, max_penalty=5)
    

def test_motif_missing_file():
    with pytest.raises(FileNotFoundError):
        MotiFind(None, "doesnt_exist.txt", max_penalty=5)
        

def test_no_match_found():
    # create object without files
    mf = MotiFind(None, None, max_penalty=5)

    # manually defined sequence
    mf.headers = ["test_seq"]
    mf.sequences = ["AAAAAAA"]

    # motif that cannot occur
    mf.motif = [
        ('char', {'C'}, 0),
        ('char', {'G'}, 0)
    ]

    # run scan
    mf.ScanSeqForMotif(chosen_idx=0)

    # check result key exists
    assert "test_seq" in mf.all_matches

    # check no-match message
    assert mf.all_matches["test_seq"] == ["NO MATCHES FOUND."]


def test_penalty_threshold():
    # create object without files
    mf = MotiFind(None, None, max_penalty=3)

    # manually defined sequence
    mf.headers = ["test_seq"]
    mf.sequences = ["ATGCCATAG"]

    # motif contains two mismatches:
    mf.motif = [
        ('char', {'A'}, 0),
        ('char', {'C'}, 1),
        ('char', {'G'}, 0),
        ('char', {'T'}, 2)
    ]

    # run scan
    mf.ScanSeqForMotif(chosen_idx=0)

    # check result key exists
    assert "test_seq" in mf.all_matches

    # expected match
    expected = (0, 4, 0.0)

    # check exactly one match was found
    assert len(mf.all_matches["test_seq"]) == 1

    # compare result
    assert mf.all_matches["test_seq"][0] == expected


def test_penalty_exceeds_threshold():
    # create object without files
    mf = MotiFind(None, None, max_penalty=2)

    # manually defined sequence
    mf.headers = ["test_seq"]
    mf.sequences = ["ATGCCATAG"]

    # motif with total mismatch penalty = 3
    mf.motif = [
        ('char', {'A'}, 0),
        ('char', {'C'}, 1),
        ('char', {'G'}, 0),
        ('char', {'T'}, 2)
    ]

    # run scan
    mf.ScanSeqForMotif(chosen_idx=0)

    # check result key exists
    assert "test_seq" in mf.all_matches

    # motif should be rejected
    assert mf.all_matches["test_seq"] == ["NO MATCHES FOUND."]


def test_partial_match():
    # create object without files
    mf = MotiFind(None, None, max_penalty=4)

    # manually defined sequence
    mf.headers = ["test_seq"]
    mf.sequences = ["ATGCCATAG"]

    # motif with one mismatch
    mf.motif = [
        ('char', {'A'}, 0),
        ('char', {'C'}, 1),
        ('char', {'G'}, 0),
        ('char', {'C'}, 0)
    ]

    # run scan
    mf.ScanSeqForMotif(chosen_idx=0)

    # expected confidence:
    expected = (0, 4, 0.75)

    # check match exists
    assert expected in mf.all_matches["test_seq"]
    