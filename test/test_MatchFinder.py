import pytest
from motifind import MotiFind

def test_scanseq_simple_match():
    # create object without files 
    mf = MotiFind(None, None, max_penalty=5)

    # manually defined FASTA-like data
    mf.headers = ["test_seq"]
    mf.sequences = ["GTCACGACTC"]

    # simple motif: A then C
    mf.motif = [
        ('char', {'A'}, 6),
        ('char', {'C'}, 6)
    ]

    # run scan
    mf.ScanSeqForMotif(chosen_idx=0)

    # check result key exists
    assert "test_seq" in mf.all_matches

    # expected matches
    expected = [
        (3, 5, 1.0),
        (6, 8, 1.0)
    ]

    # check number of matches
    assert len(mf.all_matches["test_seq"]) == 2

    # compare exact results
    assert set(mf.all_matches["test_seq"]) == set(expected)


def test_no_seq_no_motif():
    # create object without files
    mf = MotiFind(None, None, max_penalty=0)
    mf.ScanFastaForMotif()

    # should not find any matches
    assert len(mf.all_matches) == 0


def test_motif_longer_than_seq():

    # create object without files
    mf = MotiFind(None, None, max_penalty=0)

    # manually defined FASTA-like data
    mf.headers = ["test_seq"]
    mf.sequences = ["AAA"]

    # motif longer than sequence
    mf.motif = [
        ('char', {'A'}, 0),
        ('char', {'A'}, 0),
        ('char', {'A'}, 0),
        ('char', {'A'}, 0),
        ('char', {'A'}, 0)
    ]

    # should not find any matches
    mf.ScanFastaForMotif()
    assert mf.all_matches == {"test_seq": ["NO MATCHES FOUND."]}


def test_maxpen_missing():

    # create object without files and omitting max penalty
    with pytest.raises(TypeError):
            mf = MotiFind(None, None)


def test_maxpen_no_int():

    # create object without files and max penalty as int
    mf = MotiFind(None, None, max_penalty="two")

    # simple motif
    mf.motif = [
        ('char', {'A'}, 5),
        ('char', {'A'}, 5),
        ('char', {'A'}, 5)
    ]

    testseq = "AAAAA"

    with pytest.raises(TypeError):
        mf.ScanSeqForMotif(input_sequence=testseq)

            

def test_rawseq_but_no_header():

    # create object without files
    mf = MotiFind(None, None, max_penalty=0)

    # simple motif
    mf.motif = [
        ('char', {'A'}, 5),
        ('char', {'A'}, 5),
        ('char', {'A'}, 5)
    ]

    testseq = "AAAAA"

    # should not find three matches
    mf.ScanSeqForMotif(input_sequence=testseq)
    assert len(mf.all_matches) == 3


def test_empty_rawseq():

    # create object without files
    mf = MotiFind(None, None, max_penalty=0)

    # simple motif
    mf.motif = [
        ('char', {'A'}, 5),
        ('char', {'A'}, 5),
        ('char', {'A'}, 5)
    ]

    testseq = ""

    # should not find three matches
    mf.ScanSeqForMotif(input_sequence=testseq, head="test_empty_seq")
    assert mf.all_matches == {"test_empty_seq": ["NO MATCHES FOUND."]}


def test_gaps():

    # create object without files
    mf = MotiFind(None, None, max_penalty=0)

    # simple motif
    mf.motif = [
        ('char', {'A'}, 5),
        ('gap', 0, 2),
        ('char', {'A'}, 5),
        ('char', {'A'}, 5)
    ]

    testseq = "AAAAA"

    # should not find three matches
    mf.ScanSeqForMotif(input_sequence=testseq, head="test_gaps")
    assert len(mf.all_matches["test_gaps"]) == 6


def test_rawseq_simple_match():
    # create object without files
    mf = MotiFind(None, None, max_penalty=5)
    
    # raw sequence string
    rawseq = "ATTTCCAGGA"
    
    # simple motif: C then A
    mf.motif = [
        ('char', {'C'}, 6),
        ('char', {'A'}, 6)
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
        ('char', {'C'}, 6),
        ('char', {'G'}, 6)
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
    expected = [(0, 4, 0.0), (1, 5, 0.0), (2, 6, 0.33), (3, 7, 1.0), (4, 8, 0.0), (5, 9, 0.0)]

    # check exactly one match was found
    assert len(mf.all_matches["test_seq"]) == 6

    # compare result
    assert mf.all_matches["test_seq"] == expected


def test_penalty_exceeds_threshold():
    # create object without files
    mf = MotiFind(None, None, max_penalty=2)

    # manually defined sequence
    mf.headers = ["test_seq"]
    mf.sequences = ["ATGCCATAG"]

    # motif with total mismatch penalty = 3
    mf.motif = [
        ('char', {'A'}, 5),
        ('char', {'C'}, 5),
        ('char', {'G'}, 5),
        ('char', {'T'}, 5)
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
    
